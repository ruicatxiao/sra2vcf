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
    parser.add_argument('--variant_list', required=True, help='Variant list file with variants of interest')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads to use (default: 8)')
    return parser.parse_args()

def parse_variant_list(variant_list_file):
    logging.info(f"Parsing variant list {variant_list_file}...")
    variants = []
    with open(variant_list_file, 'r') as f:
        header = f.readline().strip().split()
        if header != ['CHR', 'LOCATION']:
            logging.error("Variant list file must have header 'CHR LOCATION'")
            sys.exit(1)
        for line in f:
            if not line.strip():
                continue
            cols = line.strip().split()
            if len(cols) != 2:
                logging.warning(f"Skipping invalid line in variant list: {line.strip()}")
                continue
            chr_region, location = cols[0], int(cols[1])
            variants.append({'chr_region': chr_region, 'location': location})
    return variants

def setup_directories():
    cwd = os.getcwd()
    fq_dir = os.path.join(cwd, 'fastq')
    bam_dir = os.path.join(cwd, 'bam')
    vcf_dir = os.path.join(cwd, 'vcf')
    output_dir = os.path.join(cwd, 'out')
    star_index_dir = os.path.join(cwd, 'STAR_INDEX')
    bam_cov_dir = os.path.join(cwd, 'bam_cov')

    for directory in [fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir, bam_cov_dir]:
        if not os.path.exists(directory):
            os.makedirs(directory)
    return fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir, bam_cov_dir

def check_genome_indexing(ref_genome, ref_gtf, threads, star_index_dir, variants):
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

    # 1d. Validate chromosomes in variant list
    logging.info(f"Validating chromosomes in variant list...")
    chr_set = set([variant['chr_region'] for variant in variants])
    chromosomes_in_fai = set()
    with open(fai_file, 'r') as f:
        for line in f:
            cols = line.strip().split('\t')
            chromosomes_in_fai.add(cols[0])

    missing_chromosomes = chr_set - chromosomes_in_fai
    if missing_chromosomes:
        logging.error(f"Error: Chromosomes {', '.join(missing_chromosomes)} not found in reference genome.")
        sys.exit(1)
    else:
        logging.info("All chromosomes in variant list are present in the reference genome.")

def parse_samplesheet(samplesheet_file):
    logging.info(f"Parsing samplesheet {samplesheet_file}...")
    df = pd.read_csv(samplesheet_file)
    samples = df.to_dict('records')
    return samples

def download_sra(sample, fq_dir, threads):
    sra_number = sample['sampleID']
    logging.info(f"Downloading SRA data for {sra_number} using prefetch...")
    sra_file = os.path.join(fq_dir, f"{sra_number}.sra")
    max_retries = 2  # Total number of retries
    retries = 0

    # Use prefetch to download the .sra file with retry logic
    while retries <= max_retries:
        try:
            subprocess.run([
                'prefetch', '--max-size', '100000000', sra_number, '-o', sra_file
            ], check=True)
            logging.info(f"Successfully downloaded {sra_number}.sra")
            break  # Exit the loop if prefetch is successful
        except subprocess.CalledProcessError as e:
            retries += 1
            if retries > max_retries:
                logging.error(f"prefetch failed for {sra_number} after {max_retries + 1} attempts: {e}")
                return
            else:
                logging.warning(f"prefetch attempt {retries} failed for {sra_number}. Retrying...")
                continue

    # Use fasterq-dump to convert .sra to FASTQ files
    logging.info(f"Converting .sra file to FASTQ using fasterq-dump for {sra_number}...")
    try:
        if sample['libraryType'] == 'Paired':
            subprocess.run([
                'fasterq-dump', sra_file,
                '--split-files',
                '--threads', str(threads),
                '--progress',
                '--temp', '.',
                '--details',
                '--outdir', fq_dir
            ], check=True)
            # Verify that the FASTQ files exist
            read1 = os.path.join(fq_dir, f"{sra_number}_1.fastq")
            read2 = os.path.join(fq_dir, f"{sra_number}_2.fastq")
            if not os.path.exists(read1) or not os.path.exists(read2):
                logging.error(f"FASTQ files for {sra_number} not found after fasterq-dump.")
                return
        else:
            subprocess.run([
                'fasterq-dump', sra_file,
                '--threads', str(threads),
                '--progress',
                '--temp', '.',
                '--details',
                '--outdir', fq_dir
            ], check=True)
            # Verify that the FASTQ file exists
            read1 = os.path.join(fq_dir, f"{sra_number}.fastq")
            if not os.path.exists(read1):
                logging.error(f"FASTQ file for {sra_number} not found after fasterq-dump.")
                return
    except subprocess.CalledProcessError as e:
        logging.error(f"fasterq-dump failed for {sra_number}: {e}")
        return

    # Remove the .sra file to save space
    os.remove(sra_file)

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
                '--three_prime_clip_R1', '10',
                '--three_prime_clip_R2', '10',
                raw_read1, raw_read2
            ], check=True)
            # Trim Galore outputs
            trimmed_read1 = os.path.join(fq_dir, f"{sample_id}_1_val_1.fq.gz")
            trimmed_read2 = os.path.join(fq_dir, f"{sample_id}_2_val_2.fq.gz")
            shutil.move(os.path.join(fq_dir, os.path.basename(raw_read1).replace('.fastq', '_val_1.fq.gz')), trimmed_read1)
            shutil.move(os.path.join(fq_dir, os.path.basename(raw_read2).replace('.fastq', '_val_2.fq.gz')), trimmed_read2)
            trimmed_reads['read1'] = trimmed_read1
            trimmed_reads['read2'] = trimmed_read2
        else:
            subprocess.run([
                'trim_galore', '--cores', str(threads),
                '--fastqc', '--gzip', '-o', fq_dir,
                '--three_prime_clip_R1', '10',
                raw_read1
            ], check=True)
            # Trim Galore outputs
            trimmed_read1 = os.path.join(fq_dir, f"{sample_id}_trimmed.fq.gz")
            shutil.move(os.path.join(fq_dir, os.path.basename(raw_read1).replace('.fastq', '_trimmed.fq.gz')), trimmed_read1)
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

def process_bam_file_for_sample(sample, bam_dir, bam_cov_dir, threads):
    sample_id = sample['sampleID']
    tech_type = sample['techType']
    molecule_type = sample['moleculeType']

    if tech_type == 'Illumina' and molecule_type == 'RNA':
        # For STAR alignment, the BAM file is already sorted and generated
        # We just need to index it if not already indexed
        bam_file = os.path.join(bam_dir, f"{sample_id}_Aligned.sortedByCoord.out.bam")
        bai_file = bam_file + '.bai'
        if not os.path.exists(bai_file):
            logging.info(f"Indexing {bam_file}...")
            subprocess.run([
                'samtools', 'index', '-@', str(threads), bam_file
            ], check=True)
    else:
        sam_file = os.path.join(bam_dir, f"{sample_id}.sam")
        bam_file = sam_file.replace('.sam', '.bam')
        if os.path.exists(sam_file):
            logging.info(f"Converting and sorting {sam_file} to {bam_file}...")
            subprocess.run([
                'samtools', 'sort', '-@', str(threads), '-o', bam_file, sam_file
            ], check=True)
            os.remove(sam_file)
            logging.info(f"Indexing {bam_file}...")
            subprocess.run([
                'samtools', 'index', '-@', str(threads), bam_file
            ], check=True)
        else:
            logging.error(f"SAM file {sam_file} not found for sample {sample_id}")
            return

    # Calculate coverage
    coverage_output = os.path.join(bam_cov_dir, f"{sample_id}_cov.txt")
    logging.info(f"Calculating coverage for {bam_file}...")
    subprocess.run([
        'samtools', 'coverage', '-o', coverage_output, '-d', '0', bam_file
    ], check=True)

def check_bam_exists(sample, bam_dir):
    sample_id = sample['sampleID']
    tech_type = sample['techType']
    molecule_type = sample['moleculeType']

    if tech_type == 'Illumina' and molecule_type == 'RNA':
        bam_file = os.path.join(bam_dir, f"{sample_id}_Aligned.sortedByCoord.out.bam")
    else:
        bam_file = os.path.join(bam_dir, f"{sample_id}.bam")

    bai_file = bam_file + '.bai'
    if os.path.exists(bam_file) and os.path.exists(bai_file):
        return True
    else:
        return False

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
        return

    if tech_type == 'Illumina' and sample['moleculeType'] == 'RNA':
        bam_file = os.path.join(bam_dir, f"{sample_id}_Aligned.sortedByCoord.out.bam")
    else:
        bam_file = os.path.join(bam_dir, f"{sample_id}.bam")

    if not os.path.exists(bam_file):
        logging.error(f"BAM file for sample {sample_id} not found.")
        return

    raw_vcf = os.path.join(vcf_dir, f"{sample_id}_raw.vcf")
    try:
        subprocess.run([
            'bcftools', 'mpileup', '--config', config,
            '--threads', str(threads), '-d', '5000', '--output-type', 'v',
            '--max-idepth', '5000',
            '--annotate', 'INFO/AD',
            '-o', raw_vcf,
            '-f', ref_genome, bam_file
        ], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Variant discovery failed for sample {sample_id}: {e}")
        return
    except Exception as e:
        logging.error(f"An error occurred during variant discovery for sample {sample_id}: {e}")
        return

def variant_calling(sample, vcf_dir, threads):
    sample_id = sample['sampleID']
    logging.info(f"Performing variant calling for {sample_id}...")
    raw_vcf = os.path.join(vcf_dir, f"{sample_id}_raw.vcf")
    called_vcf = os.path.join(vcf_dir, f"{sample_id}_called.vcf")
    final_vcf = os.path.join(vcf_dir, f"{sample_id}_final.vcf")

    if not os.path.exists(raw_vcf):
        logging.error(f"Raw VCF file for sample {sample_id} not found.")
        return

    try:
        subprocess.run([
            'bcftools', 'call', '--threads', str(threads),
            '-A', '-mv',
            '-o', called_vcf,
            raw_vcf
        ], check=True)
        subprocess.run([
            'bcftools', 'filter', '--output-type', 'v', '--threads', str(threads),
            '-i', 'QUAL>=20 && INFO/DP>=20',
            '--output', final_vcf,
            called_vcf
        ], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Variant calling failed for sample {sample_id}: {e}")
        return
    except Exception as e:
        logging.error(f"An error occurred during variant calling for sample {sample_id}: {e}")
        return

def analyze_results(samples, vcf_dir, output_dir, variants):
    for variant in variants:
        chr_region = variant['chr_region']
        location = variant['location']
        output_file = os.path.join(output_dir, f"{chr_region}_{location}_output.tsv")
        logging.info(f"Analyzing results for {chr_region}:{location} and writing to {output_file}...")

        with open(output_file, 'w') as out_f:
            out_f.write('sampleID\tmoleculeType\tlibraryType\ttechType\tIDV\tDP\tIMF\tAD\tAF\n')

            for sample in samples:
                sample_id = sample['sampleID']
                final_vcf = os.path.join(vcf_dir, f"{sample_id}_final.vcf")

                if not os.path.exists(final_vcf):
                    logging.warning(f"Final VCF for sample {sample_id} not found. Skipping.")
                    continue

                with open(final_vcf, 'r') as vcf_f:
                    for line in vcf_f:
                        if line.startswith('#'):
                            continue  # Skip header lines

                        cols = line.strip().split('\t')

                        # Extract necessary VCF columns
                        chrom = cols[0]
                        pos = int(cols[1])
                        ref = cols[3]
                        alt = cols[4]
                        info = cols[7]

                        # Check if the variant matches the specified variant
                        if chrom == chr_region and pos == location:
                            # Parse INFO fields into a dictionary
                            info_fields = info.split(';')
                            info_dict = {}
                            for field in info_fields:
                                if '=' in field:
                                    key, value = field.split('=', 1)
                                    info_dict[key] = value
                                else:
                                    # For flags like 'INDEL' without '=', set the key with value None
                                    info_dict[field] = None

                            # Extract common fields
                            idv = info_dict.get('IDV', 'NA')
                            dp = info_dict.get('DP', 'NA')
                            imf = info_dict.get('IMF', 'NA')

                            # Determine if the variant is an INDEL or SNP
                            if 'INDEL' in info_dict:
                                # INDEL Variant
                                ad = 'NA'
                                af = 'NA'
                            else:
                                # SNP Variant
                                # Extract AD and calculate AF
                                ad = info_dict.get('AD', 'NA')
                                if ad != 'NA':
                                    ad_values = ad.split(',')
                                    if len(ad_values) >= 2:
                                        try:
                                            # Sum all alternate allele counts (excluding the first value which is for the reference)
                                            ad_alts = sum(float(ad_value) for ad_value in ad_values[1:])
                                            dp_val = float(dp) if dp != 'NA' else 0
                                            af = ad_alts / dp_val if dp_val > 0 else 'NA'
                                            af = f"{af:.4f}" if isinstance(af, float) else 'NA'
                                            # Format AD to include all allele depths (reference and alternates)
                                            ad = ','.join(ad_values)
                                        except ValueError:
                                            ad = 'NA'
                                            af = 'NA'
                                    else:
                                        ad = 'NA'
                                        af = 'NA'
                                else:
                                    ad = 'NA'
                                    af = 'NA'

                            # Write the results to the output TSV
                            out_f.write(f"{sample_id}\t{sample['moleculeType']}\t{sample['libraryType']}\t{sample['techType']}\t{idv}\t{dp}\t{imf}\t{ad}\t{af}\n")

                            break  # Assuming only one matching variant per sample

        # Generate plots for this variant
        tsv_file = output_file
        plot_results_with_r(output_dir, tsv_file, f"{chr_region}_{location}")

def plot_results_with_r(output_dir, tsv_file, output_prefix):
    logging.info(f"Generating plots for {output_prefix} using R script...")
    bin_dir = os.path.join(os.getcwd(), 'bin')
    r_script = os.path.join(bin_dir, 'plot_variants.R')
    if not os.path.exists(r_script):
        logging.error(f"R script {r_script} not found in bin directory.")
        return
    try:
        # Pass the output directory and output prefix to the R script
        subprocess.run(['Rscript', r_script, tsv_file, output_dir, output_prefix], check=True)
        logging.info(f"Plots generated successfully for {output_prefix}.")
    except subprocess.CalledProcessError as e:
        logging.error(f"R script failed for {output_prefix}: {e}")

def run_multiqc(fq_dir, bam_cov_dir):
    multiqc_output_dir = os.path.join(os.getcwd(), 'multiqc_output')
    if not os.path.exists(multiqc_output_dir):
        os.makedirs(multiqc_output_dir)
    logging.info("Running multiQC...")
    try:
        subprocess.run([
            'multiqc', '--outdir', multiqc_output_dir, fq_dir, bam_cov_dir
        ], check=True)
        logging.info("multiQC analysis complete.")
    except subprocess.CalledProcessError as e:
        logging.error(f"multiQC failed: {e}")

def main():
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')

    args = parse_args()

    # Setup directories
    fq_dir, bam_dir, vcf_dir, output_dir, star_index_dir, bam_cov_dir = setup_directories()

    # Copy reference files to current working directory
    cwd = os.getcwd()
    ref_genome = os.path.join(cwd, 'reference_genome.fasta')
    ref_gtf = os.path.join(cwd, 'reference_genome.gtf')
    if not os.path.exists(ref_genome):
        shutil.copy(args.ref_genome, ref_genome)
    if not os.path.exists(ref_gtf):
        shutil.copy(args.gtf, ref_gtf)

    # Parse variant list
    variants = parse_variant_list(args.variant_list)

    # Check genome indexing
    check_genome_indexing(ref_genome, ref_gtf, args.threads, star_index_dir, variants)

    # Parse samplesheet
    samples = parse_samplesheet(args.sample)

    # Process each sample
    for sample in samples:
        sample_type = sample['sampleType']
        tech_type = sample['techType']
        library_type = sample['libraryType']
        logging.info(f"Processing sample {sample['sampleID']} ({sample_type}, {tech_type}, {library_type})")

        # Check if BAM and .bai files exist
        bam_exists = check_bam_exists(sample, bam_dir)

        if bam_exists:
            logging.info(f"BAM files for sample {sample['sampleID']} already exist. Skipping steps I to III.")
        else:
            # Steps I to III
            trimmed_reads = {}
            # Step II: Branching
            if sample_type == 'SRA':
                download_sra(sample, fq_dir, args.threads)
                trimmed_reads = trim_reads(sample, fq_dir, args.threads)
            elif sample_type == 'LOCAL':
                trimmed_reads = trim_reads(sample, fq_dir, args.threads)
            # Alignment
            align_reads(sample, fq_dir, bam_dir, ref_genome, star_index_dir, args.threads, trimmed_reads)
            # BAM processing
            process_bam_file_for_sample(sample, bam_dir, bam_cov_dir, args.threads)

    # Step IV: Variant discovery and calling
    for sample in samples:
        sample_id = sample['sampleID']
        final_vcf = os.path.join(vcf_dir, f"{sample_id}_final.vcf")
        if os.path.exists(final_vcf) and os.path.getsize(final_vcf) > 0:
            logging.info(f"Final VCF for sample {sample_id} already exists and is not empty. Skipping variant discovery and calling.")
            continue
        variant_discovery(sample, bam_dir, vcf_dir, ref_genome, args.threads)
        variant_calling(sample, vcf_dir, args.threads)

    # Step V: Analysis and plotting for each variant
    analyze_results(samples, vcf_dir, output_dir, variants)

    # Run multiQC
    run_multiqc(fq_dir, bam_cov_dir)

if __name__ == '__main__':
    main()
