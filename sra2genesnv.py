import os
import sys
import subprocess
import pandas as pd
import argparse
import shutil
import logging

# Argument Parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Variant Calling Pipeline")
    parser.add_argument("--ref_genome", required=True, help="Path to reference genome (.fa/.fasta/.fna)")
    parser.add_argument("--gtf", required=True, help="Path to reference GTF file")
    parser.add_argument("--sample", required=True, help="Path to samplesheet.csv")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use (default: 8)")
    parser.add_argument("--chr_region", required=True, help="Chromosome region of interest")
    parser.add_argument("--region_start", type=int, required=True, help="Region start position")
    parser.add_argument("--region_end", type=int, required=True, help="Region end position")
    return parser.parse_args()

# Utility Functions
def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"Error executing command: {e}")

# Check and create directories
def create_directories():
    for directory in ["fastq", "bam", "vcf", "out", "STAR_INDEX"]:
        if not os.path.exists(directory):
            os.makedirs(directory)

# Indexing Functions
def index_reference_genome(ref_genome, ref_gtf, threads):
    # 1a. Samtools faidx
    if not os.path.exists(f"{ref_genome}.fai"):
        run_command(f"samtools faidx {ref_genome}")

    # 1b. BWA index
    bwa_extensions = [".sa", ".amb", ".ann", ".pac", ".bwt"]
    if not all(os.path.exists(ref_genome + ext) for ext in bwa_extensions):
        run_command(f"bwa index {ref_genome}")

    # 1c. STAR index
    if not os.listdir("STAR_INDEX"):
        run_command(
            f"STAR --runThreadN {threads} --runMode genomeGenerate --sjdbGTFfile {ref_gtf} \
            --genomeFastaFiles {ref_genome} --genomeDir STAR_INDEX --genomeSAindexNbases 10"
        )

# Validate region boundaries
def validate_region(ref_genome, chr_region, region_start, region_end):
    fai_file = f"{ref_genome}.fai"
    with open(fai_file, 'r') as f:
        for line in f:
            if line.startswith(chr_region):
                length = int(line.split()[1])
                if 1 <= region_start < region_end <= length:
                    return True
                else:
                    sys.exit("Error: Region boundaries are invalid.")
    sys.exit("Error: Chromosome region not found in reference genome.")

# SRA Download Function
def download_sra(sample_id, threads):
    run_command(f"fasterq-dump -e {threads} --split-files -p {sample_id} -O fastq")

# Trimming Function
def trim_reads(read1, read2, threads, paired=True):
    if paired:
        run_command(f"trim_galore --cores {threads} --paired --fastqc --gzip -o fastq {read1} {read2}")
    else:
        run_command(f"trim_galore --cores {threads} --fastqc --gzip -o fastq {read1}")

# Alignment Function
def align_reads(ref_genome, read1, read2, sample_name, library_type, tech_type, threads):
    bam_dir = "bam"
    if tech_type == "Illumina" and library_type == "Paired":
        run_command(f"bwa mem -t {threads} {ref_genome} {read1} {read2} > {bam_dir}/{sample_name}.sam")
    elif tech_type == "Illumina" and library_type == "Single":
        run_command(f"bwa mem -t {threads} {ref_genome} {read1} > {bam_dir}/{sample_name}.sam")
    elif tech_type in ["Pacbio", "ONT"]:
        mode = "map-pb" if tech_type == "Pacbio" else "map-ont"
        run_command(f"minimap2 -ax {mode} -t {threads} --secondary=no --sam-hit-only {ref_genome} {read1} > {bam_dir}/{sample_name}.sam")
    elif tech_type == "Illumina" and library_type == "RNA":
        bam_prefix = f"{bam_dir}/{sample_name}_"
        run_command(
            f"STAR --runMode alignReads --runThreadN {threads} --genomeDir STAR_INDEX \
            --readFilesIn {read1} {read2 if read2 else ''} --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts --bamRemoveDuplicatesType UniqueIdentical --twopassMode Basic \
            --outFileNamePrefix {bam_prefix} --outWigStrand Stranded --outSAMstrandField intronMotif \
            --alignIntronMax 2500 --chimOutType Junctions --limitBAMsortRAM 1000000000 --readFilesCommand gunzip -c"
        )

# BAM Processing Function
def process_bam_files(threads):
    bam_dir = "bam"
    for sam_file in os.listdir(bam_dir):
        if sam_file.endswith(".sam"):
            bam_file = sam_file.replace(".sam", ".bam")
            run_command(f"samtools view -h -@ {threads} -S -b {bam_dir}/{sam_file} | samtools sort -@ {threads} -o {bam_dir}/{bam_file}")
            os.remove(os.path.join(bam_dir, sam_file))
            run_command(f"samtools index -@ {threads} {bam_dir}/{bam_file}")

# Variant Calling Function
def variant_calling(ref_genome, bam_dir, vcf_dir, tech_type, threads):
    for bam_file in os.listdir(bam_dir):
        if bam_file.endswith(".bam"):
            sample_name = bam_file.replace(".bam", "")
            vcf_raw = f"{vcf_dir}/{sample_name}_raw.vcf"
            config = "illumina-1.20" if tech_type == "Illumina" else "pacbio-ccs-1.20" if tech_type == "Pacbio" else "ont-sup-1.20"
            run_command(f"bcftools mpileup --config {config} --threads {threads} -d 5000 --output-type v --max-idepth 5000 -o {vcf_raw} -f {ref_genome} {bam_dir}/{bam_file}")
            vcf_called = vcf_raw.replace("_raw.vcf", "_called.vcf")
            run_command(f"bcftools call --threads {threads} -A -mv -o {vcf_called} {vcf_raw}")
            vcf_final = vcf_called.replace("_called.vcf", "_final.vcf")
            run_command(f"bcftools filter --output-type v --threads {threads} --output {vcf_final} {vcf_called}")

# Analyze Variants and Output Results
def analyze_variants(chr_region, region_start, region_end, vcf_dir, output_dir, samplesheet):
    output_file = f"{output_dir}/{chr_region}_{region_start}_{region_end}_output.tsv"
    with open(output_file, 'w') as out_f:
        out_f.write("sampleID,moleculeType,libraryType,techType,IDV,DP,INFper\n")
        for vcf_file in os.listdir(vcf_dir):
            if vcf_file.endswith("_final.vcf"):
                with open(f"{vcf_dir}/{vcf_file}", 'r') as vcf_f:
                    for line in vcf_f:
                        if not line.startswith("#"):
                            cols = line.strip().split('\t')
                            if cols[0] == chr_region and region_start <= int(cols[1]) <= region_end:
                                info_fields = cols[7].split(';')
                                info_dict = {field.split('=')[0]: field.split('=')[1] for field in info_fields if '=' in field}
                                idv = info_dict.get('IDV', 'NA')
                                dp = info_dict.get('DP', 'NA')
                                inf = info_dict.get('IMF', 'NA')
                                sample_id = vcf_file.split("_final")[0]
                                sample_details = samplesheet.loc[samplesheet['sampleID'] == sample_id].iloc[0]
                                molecule_type = sample_details['moleculeType']
                                library_type = sample_details['libraryType']
                                tech_type = sample_details['techType']
                                out_f.write(f"{sample_id},{molecule_type},{library_type},{tech_type},{idv},{dp},{inf}\n")

# Plotting Function
def plot_indel_frequency(output_dir, chr_region, region_start, region_end):
    import matplotlib.pyplot as plt
    import pandas as pd

    tsv_file = f"{output_dir}/{chr_region}_{region_start}_{region_end}_output.tsv"
    data = pd.read_csv(tsv_file, sep=',')
    plt.figure(figsize=(10, 6))
    for idx, row in data.iterrows():
        idv = float(row['IDV']) if row['IDV'] != 'NA' else 0
        dp = float(row['DP']) if row['DP'] != 'NA' else 0
        deletion = dp - idv
        plt.bar(row['sampleID'], [idv, deletion], label=['Reference', 'Deletion'], stacked=True)

    plt.title(f"INDEL Frequency in Region {chr_region}:{region_start}-{region_end}")
    plt.xlabel("Sample ID")
    plt.ylabel("Frequency")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{chr_region}_{region_start}_{region_end}_plot.png")
    plt.close()

# Main Function
def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s:%(message)s')
    args = parse_args()
    create_directories()

    samplesheet = pd.read_csv(args.sample)
    ref_genome = args.ref_genome
    ref_gtf = args.gtf
    threads = args.threads
    chr_region = args.chr_region
    region_start = args.region_start
    region_end = args.region_end

    index_reference_genome(ref_genome, ref_gtf, threads)
    validate_region(ref_genome, chr_region, region_start, region_end)

    for _, sample in samplesheet.iterrows():
        if sample['sampleType'] == "SRA":
            download_sra(sample['sampleID'], threads)
        # Assume reads are trimmed and available in 'fastq'
        read1, read2 = sample['read1'], sample['read2']
        align_reads(ref_genome, read1, read2, sample['sampleID'], sample['libraryType'], sample['techType'], threads)

    process_bam_files(threads)
    variant_calling(ref_genome, "bam", "vcf", "Illumina", threads)
    analyze_variants(chr_region, region_start, region_end, "vcf", "out", samplesheet)
    plot_indel_frequency("out", chr_region, region_start, region_end)

if __name__ == "__main__":
    main()