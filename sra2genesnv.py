#!/usr/bin/env python3

import argparse
import os
import sys
import subprocess
import pandas as pd
import shutil
import logging

def parse_args():
    parser = argparse.ArgumentParser(description="Variant Calling Pipeline")
    parser.add_argument('--ref_genome', required=True, help='Reference genome file (FASTA format)')
    parser.add_argument('--gtf', required=True, help='Reference annotation GTF file')
    parser.add_argument('--sample', required=True, help='Samplesheet CSV file')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use (default: 8)')
    parser.add_argument('--chr_region', required=True, help='Chromosome region of interest')
    parser.add_argument('--region_start', type=int, required=True, help='Region start position')
    parser.add_argument('--region_end', type=int, required=True, help='Region end position')
    return parser.parse_args()

def setup_directories():
    cwd = os.getcwd()
    fq_dir = os.path.join(cwd, 'fastq')
    bam_dir = os.path.join(cwd, 'bam')
    vcf_dir = os.path.join(cwd, 'vcf')
    output_dir = os.path.join(cwd, 'out')
    star_index_dir = os.path.join(cwd, 'STAR_INDEX')

    for directory in [fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir]:
        if not os.path.exists(directory):
            os.makedirs(directory)
    return fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir

def check_genome_indexing(ref_genome, ref_gtf, threads, star_index_dir, chr_region, region_start, region_end):
    # 1a. Check if genome has been faidx processed
    fai_file = ref_genome + '.fai'
    if not os.path.exists(fai_file):
        logging.info(f"Indexing reference genome with samtools faidx...")
        subprocess.run(['samtools', 'faidx', ref_genome], check=True)
    else:
        logging.info(f"FAI index file {fai_file} already exists. Skipping samtools faidx.")

    # 1b. Check if genome has been BWA indexed
    bwa_index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    bwa_index_files = [ref_genome + ext for ext in bwa_index_extensions]
    if not all(os.path.exists(f) for f in bwa_index_files):
        logging.info(f"Indexing reference genome with bwa index...")
        subprocess.run(['bwa', 'index', ref_genome], check=True)
    else:
        logging.info(f"BWA index files already exist. Skipping bwa index.")

    # 1c. Check if genome has been STAR indexed
    if not os.listdir(star_index_dir):
        logging.info(f"Generating STAR genome index...")
        subprocess.run([
            'STAR', '--runThreadN', str(threads),
            '--runMode', 'genomeGenerate',
            '--sjdbGTFfile', ref_gtf,
            '--genomeFastaFiles', ref_genome,
            '--genomeDir', star_index_dir,
            '--genomeSAindexNbases', '10'
        ], check=True)
    else:
        logging.info(f"STAR genome index already exists in {star_index_dir}. Skipping STAR genomeGenerate.")

    # 1d. Validate chromosome region
    logging.info(f"Validating chromosome region {chr_region}:{region_start}-{region_end}...")
    chr_valid = False
    chr_length = None
    with open(fai_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            if cols[0] == chr_region:
                chr_valid = True
                chr_length = int(cols[1])
                break
    if not chr_valid:
        logging.error(f"Error: Chromosome {chr_region} not found in reference genome.")
        sys.exit(1)
    if not (1 <= region_start < region_end <= chr_length):
        logging.error(f"Error: Invalid region coordinates. Must satisfy 1 <= region_start < region_end <= {chr_length}.")
        sys.exit(1)
    logging.info(f"Chromosome region is valid.")

def parse_samplesheet(samplesheet_file):
    logging.info(f"Parsing samplesheet {samplesheet_file}...")
    df = pd.read_csv(samplesheet_file)
    samples = df.to_dict('records')
    return samples

def download_sra(sample, fq_dir, threads):
    sra_number = sample['sampleID']
    logging.info(f"Downloading SRA data for {sra_number}...")
    if sample['libraryType'] == 'Paired':
        subprocess.run([
            'fasterq-dump', '-e', str(threads),
            '--split-files', '-p', sra_number,
            '-O', fq_dir
        ], check=True)
    else:
        subprocess.run([
            'fasterq-dump', '-e', str(threads),
            '-p', sra_number,
            '-O', fq_dir
        ], check=True)

def trim_reads(sample, fq_dir, threads):
    sample_id = sample['sampleID']
    sample_type = sample['sampleType']
    library_type = sample['libraryType']
    tech_type = sample['techType']
    logging.info(f"Trimming reads for {sample_id}...")
    trimmed_reads = {}
    if sample_type == 'SRA':
        if library_type == 'Paired':
            raw_read1 = os.path.join(fq_dir, f"{sample_id}_1.fastq")
            raw_read2 = os.path.join(fq_dir, f"{sample_id}_2.fastq")
        else:
            raw_read1 = os.path.join(fq_dir, f"{sample_id}.fastq")
            raw_read2 = None
    else:  # LOCAL
        # Copy local reads to fq_dir if not already there
        if library_type == 'Paired':
            raw_read1 = os.path.join(fq_dir, os.path.basename(sample['read1']))
            raw_read2 = os.path.join(fq_dir, os.path.basename(sample['read2']))
            if not os.path.exists(raw_read1):
                shutil.copy(sample['read1'], raw_read1)
            if not os.path.exists(raw_read2):
                shutil.copy(sample['read2'], raw_read2)
        else:
            raw_read1 = os.path.join(fq_dir, os.path.basename(sample['read1']))
            raw_read2 = None
            if not os.path.exists(raw_read1):
                shutil.copy(sample['read1'], raw_read1)

    # For techType == 'Illumina', perform trimming
    if tech_type == 'Illumina':
        if library_type == 'Paired':
            subprocess.run([
                'trim_galore', '--cores', str(threads),
                '--paired', '--fastqc', '--gzip', '-o', fq_dir,
                raw_read1, raw_read2
            ], check=True)
            # Trim Galore outputs
            trimmed_read1 = os.path.join(fq_dir, os.path.basename(raw_read1).replace('.fastq', '_val_1.fq.gz'))
            trimmed_read2 = os.path.join(fq_dir, os.path.basename(raw_read2).replace('.fastq', '_val_2.fq.gz'))
            trimmed_reads['read1'] = trimmed_read1
            trimmed_reads['read2'] = trimmed_read2
        else:
            subprocess.run([
                'trim_galore', '--cores', str(threads),
                '--fastqc', '--gzip', '-o', fq_dir,
                raw_read1
            ], check=True)
            # Trim Galore outputs
            trimmed_read1 = os.path.join(fq_dir, os.path.basename(raw_read1).replace('.fastq', '_trimmed.fq.gz'))
            trimmed_reads['read1'] = trimmed_read1
    else:
        # For non-Illumina, no trimming is performed
        if library_type == 'Paired':
            trimmed_reads['read1'] = raw_read1
            trimmed_reads['read2'] = raw_read2
        else:
            trimmed_reads['read1'] = raw_read1
    return trimmed_reads

def align_reads(sample, fq_dir, bam_dir, ref_genome, star_index_dir, threads, trimmed_reads):
    sample_id = sample['sampleID']
    molecule_type = sample['moleculeType']
    library_type = sample['libraryType']
    tech_type = sample['techType']
    logging.info(f"Aligning reads for {sample_id}...")

    # Determine input reads
    read1 = trimmed_reads.get('read1')
    read2 = trimmed_reads.get('read2')

    # Alignment commands
    if tech_type == 'Illumina' and molecule_type == 'DNA':
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        if library_type == 'Paired':
            subprocess.run([
                'bwa', 'mem', '-t', str(threads),
                ref_genome, read1, read2
            ], stdout=open(sam_file, 'w'), check=True)
        else:
            subprocess.run([
                'bwa', 'mem', '-t', str(threads),
                ref_genome, read1
            ], stdout=open(sam_file, 'w'), check=True)
    elif tech_type == 'Illumina' and molecule_type == 'RNA':
        bam_prefix = os.path.join(bam_dir, f"{sample_id}_")
        if library_type == 'Paired':
            subprocess.run([
                'STAR', '--runMode', 'alignReads', '--runThreadN', str(threads),
                '--genomeDir', star_index_dir,
                '--readFilesIn', read1, read2,
                '--outSAMtype', 'BAM', 'SortedByCoordinate', '--quantMode', 'GeneCounts',
                '--bamRemoveDuplicatesType', 'UniqueIdentical', '--twopassMode', 'Basic',
                '--outFileNamePrefix', bam_prefix,
                '--outWigStrand', 'Stranded', '--outSAMstrandField', 'intronMotif',
                '--alignIntronMax', '2500', '--chimOutType', 'Junctions',
                '--limitBAMsortRAM', '1000000000', '--readFilesCommand', 'gunzip', '-c'
            ], check=True)
        else:
            subprocess.run([
                'STAR', '--runMode', 'alignReads', '--runThreadN', str(threads),
                '--genomeDir', star_index_dir,
                '--readFilesIn', read1,
                '--outSAMtype', 'BAM', 'SortedByCoordinate', '--quantMode', 'GeneCounts',
                '--bamRemoveDuplicatesType', 'UniqueIdentical', '--twopassMode', 'Basic',
                '--outFileNamePrefix', bam_prefix,
                '--outWigStrand', 'Stranded', '--outSAMstrandField', 'intronMotif',
                '--alignIntronMax', '2500', '--chimOutType', 'Junctions',
                '--limitBAMsortRAM', '1000000000', '--readFilesCommand', 'gunzip', '-c'
            ], check=True)
    elif tech_type == 'Pacbio' and molecule_type == 'DNA':
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        subprocess.run([
            'minimap2', '-ax', 'map-pb', '-t', str(threads),
            '--secondary=no', '--sam-hit-only',
            ref_genome, read1
        ], stdout=open(sam_file, 'w'), check=True)
    elif tech_type == 'ONT' and molecule_type == 'DNA':
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        subprocess.run([
            'minimap2', '-ax', 'map-ont', '-t', str(threads),
            '--secondary=no', '--sam-hit-only',
            ref_genome, read1
        ], stdout=open(sam_file, 'w'), check=True)
    elif tech_type == 'Pacbio' and molecule_type == 'RNA':
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        subprocess.run([
            'minimap2', '-ax', 'splice:hq', '-t', str(threads),
            '-uf', '-G', '2500', '--secondary=no', '--sam-hit-only',
            ref_genome, read1
        ], stdout=open(sam_file, 'w'), check=True)
    elif tech_type == 'ONT' and molecule_type == 'RNA':
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        subprocess.run([
            'minimap2', '-ax', 'splice', '-t', str(threads),
            '-uf', '-G', '2500', '-k14', '--secondary=no', '--sam-hit-only',
            ref_genome, read1
        ], stdout=open(sam_file, 'w'), check=True)
    else:
        logging.error(f"Unsupported combination for sample {sample_id}.")
        sys.exit(1)

def process_bam_files(bam_dir, threads):
    logging.info("Processing SAM to BAM and indexing...")
    for filename in os.listdir(bam_dir):
        if filename.endswith('.sam'):
            sam_path = os.path.join(bam_dir, filename)
            bam_path = sam_path.replace('.sam', '.bam')
            logging.info(f"Converting and sorting {sam_path} to {bam_path}...")
            subprocess.run([
                'samtools', 'sort', '-@', str(threads), '-o', bam_path, sam_path
            ], check=True)
            os.remove(sam_path)
    # Index BAM files
    for filename in os.listdir(bam_dir):
        if filename.endswith('.bam'):
            bam_path = os.path.join(bam_dir, filename)
            logging.info(f"Indexing {bam_path}...")
            subprocess.run([
                'samtools', 'index', '-@', str(threads), bam_path
            ], check=True)

def variant_discovery(sample, bam_dir, vcf_dir, ref_genome, threads):
    sample_id = sample['sampleID']
    tech_type = sample['techType']
    logging.info(f"Performing variant discovery for {sample_id}...")
    if tech_type == 'Illumina':
        config = 'illumina-1.20'
    elif tech_type == 'Pacbio':
        config = 'pacbio-ccs-1.20'
    elif tech_type == 'ONT':
        config = 'ont-sup-1.20'
    else:
        logging.error(f"Unsupported tech type {tech_type} for sample {sample_id}.")
        sys.exit(1)
    bam_file = None
    # For RNA samples aligned with STAR, BAM files are named differently
    if tech_type == 'Illumina' and sample['moleculeType'] == 'RNA':
        bam_file = os.path.join(bam_dir, f"{sample_id}_Aligned.sortedByCoord.out.bam")
    else:
        bam_file = os.path.join(bam_dir, f"{sample_id}.bam")
    if not os.path.exists(bam_file):
        logging.error(f"BAM file for sample {sample_id} not found.")
        sys.exit(1)
    raw_vcf = os.path.join(vcf_dir, f"{sample_id}_raw.vcf")
    subprocess.run([
        'bcftools', 'mpileup', '--config', config,
        '--threads', str(threads), '-d', '5000', '--output-type', 'v',
        '--max-idepth', '5000',
        '-o', raw_vcf,
        '-f', ref_genome, bam_file
    ], check=True)

def variant_calling(sample, vcf_dir, threads):
    sample_id = sample['sampleID']
    logging.info(f"Performing variant calling for {sample_id}...")
    raw_vcf = os.path.join(vcf_dir, f"{sample_id}_raw.vcf")
    called_vcf = os.path.join(vcf_dir, f"{sample_id}_called.vcf")
    final_vcf = os.path.join(vcf_dir, f"{sample_id}_final.vcf")
    subprocess.run([
        'bcftools', 'call', '--threads', str(threads),
        '-A', '-mv',
        '-o', called_vcf,
        raw_vcf
    ], check=True)
    subprocess.run([
        'bcftools', 'filter', '--output-type', 'v', '--threads', str(threads),
        '--output', final_vcf,
        called_vcf
    ], check=True)

def analyze_results(samples, vcf_dir, output_dir, chr_region, region_start, region_end):
    output_file = os.path.join(output_dir, f"{chr_region}_{region_start}_{region_end}_output.tsv")
    logging.info(f"Analyzing results and writing to {output_file}...")
    with open(output_file, 'w') as out_f:
        out_f.write('sampleID\tmoleculeType\tlibraryType\ttechType\tIDV\tDP\tIMF\n')
        for sample in samples:
            sample_id = sample['sampleID']
            final_vcf = os.path.join(vcf_dir, f"{sample_id}_final.vcf")
            if not os.path.exists(final_vcf):
                logging.warning(f"Final VCF for sample {sample_id} not found. Skipping.")
                continue
            with open(final_vcf, 'r') as vcf_f:
                for line in vcf_f:
                    if line.startswith('#'):
                        continue
                    cols = line.strip().split('\t')
                    chrom = cols[0]
                    pos = int(cols[1])
                    if chrom == chr_region and region_start <= pos <= region_end:
                        info_fields = cols[7].split(';')
                        info_dict = {}
                        for field in info_fields:
                            if '=' in field:
                                key, value = field.split('=')
                                info_dict[key] = value
                        idv = info_dict.get('IDV', 'NA')
                        dp = info_dict.get('DP', 'NA')
                        imf = info_dict.get('IMF', 'NA')
                        out_f.write(f"{sample_id}\t{sample['moleculeType']}\t{sample['libraryType']}\t{sample['techType']}\t{idv}\t{dp}\t{imf}\n")
                        break  # Assuming only one variant per sample in region

def plot_results(output_dir, chr_region, region_start, region_end):
    import matplotlib.pyplot as plt
    import pandas as pd
    output_file = os.path.join(output_dir, f"{chr_region}_{region_start}_{region_end}_output.tsv")
    df = pd.read_csv(output_file, sep='\t')
    if df.empty:
        logging.warning("No data to plot.")
        return
    samples = df['sampleID']
    idv = df['IDV'].astype(float)
    dp = df['DP'].astype(float)
    deletion = dp - idv
    data = pd.DataFrame({'Reference': idv, 'Deletion': deletion}, index=samples)
    data_pct = data.div(data.sum(axis=1), axis=0)
    ax = data_pct.plot(kind='bar', stacked=True)
    plt.ylabel('Percentage')
    plt.title('INDEL Frequency')
    plt.legend(loc='upper right')
    plt.tight_layout()
    plot_file = os.path.join(output_dir, f"{chr_region}_{region_start}_{region_end}_plot.png")
    plt.savefig(plot_file)
    logging.info(f"Plot saved to {plot_file}")

def main():
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

    args = parse_args()

    # Setup directories
    fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir = setup_directories()

    # Copy reference files to current working directory
    cwd = os.getcwd()
    ref_genome = os.path.join(cwd, 'reference_genome.fasta')
    ref_gtf = os.path.join(cwd, 'reference_genome.gtf')
    shutil.copy(args.ref_genome, ref_genome)
    shutil.copy(args.gtf, ref_gtf)

    # Check genome indexing
    check_genome_indexing(ref_genome, ref_gtf, args.threads, star_index_dir, args.chr_region, args.region_start, args.region_end)

    # Parse samplesheet
    samples = parse_samplesheet(args.sample)

    # Process each sample
    for sample in samples:
        sample_type = sample['sampleType']
        tech_type = sample['techType']
        library_type = sample['libraryType']
        logging.info(f"Processing sample {sample['sampleID']} ({sample_type}, {tech_type}, {library_type})")

        trimmed_reads = {}
        # Step II: Branching
        if sample_type == 'SRA':
            download_sra(sample, fq_dir, args.threads)
            trimmed_reads = trim_reads(sample, fq_dir, args.threads)
        elif sample_type == 'LOCAL':
            trimmed_reads = trim_reads(sample, fq_dir, args.threads)
        # Alignment
        align_reads(sample, fq_dir, bam_dir, ref_genome, star_index_dir, args.threads, trimmed_reads)

    # Step III: Shared BAM processing
    process_bam_files(bam_dir, args.threads)

    # Step IV: Variant discovery and calling
    for sample in samples:
        variant_discovery(sample, bam_dir, vcf_dir, ref_genome, args.threads)
        variant_calling(sample, vcf_dir, args.threads)

    # Step V: Analysis and plotting
    analyze_results(samples, vcf_dir, output_dir, args.chr_region, args.region_start, args.region_end)
    plot_results(output_dir, args.chr_region, args.region_start, args.region_end)

if __name__ == '__main__':
    main()
