![Image](https://github.com/user-attachments/assets/ef863606-31ec-4971-ba10-a29a4e3639a0)

# sra2vcf: A Comprehensive Variant Discovery Pipeline

sra2vcf is a Python-based pipeline for performing robust and comprehensive SNP/Indel analysis across diverse sequencing data types specifically tailed for Cryptosporidium parvum. This versatile tool can process both local files and SRA data from NCBI, supporting short and long reads, DNA and RNA sequencing data.

## Features

- **Versatile Input Support**:
  - Process data directly from NCBI's SRA database or local files
  - Support for multiple sequencing technologies (Illumina, PacBio, Oxford Nanopore)
  - Compatible with both DNA and RNA sequencing data
  - Handles paired-end, single-end, and long-read libraries

- **Efficient Processing**:
  - Skips redundant processing for previously analyzed samples
  - Multi-threaded execution for faster analysis
  - Appropriate aligner selection based on data type

- **Comprehensive Variant Analysis**:
  - Identifies both SNPs and INDELs
  - Focuses analysis on specific genomic positions of interest
  - Generates detailed metrics including depth, allele frequency, and more

- **Quality Control and Visualization**:
  - Automatic quality control through MultiQC
  - Custom plots for variant frequency comparison across samples
  - Technology-specific variant calling parameters

- **Accumulative processing 
  - Adding new samples over time without reprocessing existing data
  - Automatically regenerates combined figures and tables when new samples are added

## Requirements

### Software Dependencies

All tools should be available in your `$PATH`:

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.10.11+ | Pipeline execution |
| SRA Toolkit | 3.0.5+ | SRA data download and extraction |
| Trim Galore | 0.6.10+ | Read quality trimming |
| BWA | 0.7.17+ | Short-read DNA alignment |
| Minimap2 | 2.26+ | Long-read alignment |
| STAR | 2.7.10b+ | RNA-seq alignment |
| samtools | 1.19.2+ | SAM/BAM manipulation |
| bcftools | 1.20+ | Variant calling and filtering |
| R | 4.4.2+ | Result visualization |
| MultiQC | 1.25.2+ | Quality control reporting |

### R Packages Required

- ggplot2
- reshape2
- dplyr
- tools

## Installation

1. Ensure all required software is installed and available in your `$PATH`
2. Clone this repository:
   ```bash
   git clone https://github.com/yourusername/sra2genesnv.git
   cd sra2genesnv
   ```
3. Make the script executable:
   ```bash
   chmod +x sra2genesnv.py
   ```

## Input Files

### 1. Reference Genome
FASTA format reference genome file. T2T Cparvum BGF annotation is provded here.

### 2. Reference GTF Annotation
Gene annotation file in GTF3 format, we recommend you clean your GTF file with latest AGAT. T2T Cparvum BGF annotation is provded here.

### 3. Variant List
Tab-delimited file containing your variants of interest:

```
CHR	LOCATION
CP141119.1  23798
CP141120.1  549739


```

### 4. Sample Sheet (CSV)
CSV file describing samples to be processed. All possible sample configurations:

```
sampleType,moleculeType,libraryType,techType,sampleID,read1,read2
SRA,DNA,Single,Illumina,SRR12345678
SRA,DNA,Paired,Illumina,SRR12345679
SRA,DNA,Long,Pacbio,SRR12345680
SRA,DNA,Long,ONT,SRR12345681
SRA,RNA,Single,Illumina,SRR12345682
SRA,RNA,Paired,Illumina,SRR12345683
SRA,RNA,Long,Pacbio,SRR12345684
SRA,RNA,Long,ONT,SRR12345685
LOCAL,DNA,Single,Illumina,sample1,/path/to/read.fastq
LOCAL,DNA,Paired,Illumina,sample2,/path/to/read1.fastq,/path/to/read2.fastq
LOCAL,DNA,Long,Pacbio,sample3,/path/to/pacbio.fastq
LOCAL,DNA,Long,ONT,sample4,/path/to/nanopore.fastq
LOCAL,RNA,Single,Illumina,sample5,/path/to/rna.fastq
LOCAL,RNA,Paired,Illumina,sample6,/path/to/rna1.fastq,/path/to/rna2.fastq
LOCAL,RNA,Long,Pacbio,sample7,/path/to/rna_pacbio.fastq
LOCAL,RNA,Long,ONT,sample8,/path/to/rna_nanopore.fastq
```

* Both relative PATH and abosolute PATH are accepatable here. 

## Usage

Basic command:

```bash
python sra2genesnv.py \
    --ref_genome reference.fasta \
    --gtf annotation.gtf \
    --sample samplesheet.csv \
    --variant_list variants.txt \
    --threads 16
```

### Arguments

| Argument | Description | Required |
|----------|-------------|----------|
| --ref_genome | Reference genome FASTA file | Yes |
| --gtf | Reference annotation GTF file | Yes |
| --sample | Sample sheet CSV file | Yes |
| --variant_list | List of variants to analyze | Yes |
| --threads | Number of threads to use (default: 8) | No |

## Pipeline Workflow

1. **Setup & Directory Creation**
   - Creates necessary directories for intermediate and final outputs
   - Copies reference files to working directory

2. **Genome Indexing**
   - Creates indices for BWA, STAR, and samtools
   - Validates chromosomes in variant list against reference

3. **Sample Processing** (for each sample)
   - **SRA Data**: Downloads and extracts FASTQ files
   - **Quality Trimming**: Applies Trim Galore for short-read data
   - **Alignment**: Uses appropriate aligner based on data type
     - BWA: short-read DNA
     - STAR: short-read RNA
     - Minimap2: PacBio/ONT (DNA/RNA)
   - **BAM Processing**: Sorting, indexing, and coverage calculation

4. **Variant Discovery & Calling**
   - BCFtools mpileup with technology-specific parameters
      - Short-reads with illumina-1.20 preset
      - Pacbio with pacbio-ccs-1.20 preset
      - ONT with ont-sup-1.20 preset
   - Variant calling 
   - Variant quality filtering (QUAL≥20, DP≥20)

5. **Result Analysis**
   - Compiles metrics for variants of interest
   - Extracts depth, frequency, and other statistics

6. **Visualization & QC**
   - Generates plots for each variant
   - Runs MultiQC for comprehensive QC report

## Accumulative Sample Processing

sra2vcf is designed to work efficiently with growing datasets:

1. **Add New Samples Anytime**
   - Simply update your sample sheet with new samples
   - Run the pipeline with the same command as before
   - Previously processed samples will be detected and skipped

2. **Smart Reprocessing**
   - Only new samples undergo the full analysis pipeline
   - Existing sample results are preserved and reused

3. **Automatic Result Aggregation**
   - All output tables are regenerated combining both previous and new samples
   - Visualization plots are updated to include all samples in the dataset
   - MultiQC reports incorporate data from the complete sample set

This feature is particularly valuable for ongoing studies where data is collected over time or for collaborative projects where multiple datasets need to be analyzed consistently.

## Directory Structure

### Input Structure

```
.
├── sra2genesnv.py
├── plot_variants.R
├── reference.fasta
├── annotation.gtf
├── variants.txt
└── samplesheet.csv
```

### Output Structure

```
.
├── STAR_INDEX/            # STAR genome indices
├── fastq/                 # Downloaded/copied FASTQ files
├── bam/                   # Alignment files (BAM)
├── bam_cov/               # Coverage statistics
├── vcf/                   # Variant calling results
│   ├── sample1_raw.vcf    # Raw variant calls
│   ├── sample1_called.vcf # Initial variant calls
│   └── sample1_final.vcf  # Filtered variants
├── out/                   # Analysis results
│   ├── CHR_POS_output.tsv # Tabulated results
│   ├── CHR_POS_plot.pdf   # Result visualization
│   └── ...
└── multiqc_output/        # QC reports
    └── multiqc_report.html
```

## Example Results

For each variant position in your variant list, the pipeline generates:

1. A TSV file with metrics for each sample
2. Visualization plots comparing variant frequencies
3. Coverage information across samples

## Troubleshooting

### Common Issues

1. **SRA Download Failures**
   - The pipeline includes retry logic for SRA downloads
   - Check your internet connection or SRA access permissions

2. **Missing Tools**
   - Ensure all required tools are properly installed and in your PATH
   - Check version compatibility if errors occur

3. **Resource Limitations**
   - Adjust the --threads parameter for your system's capabilities
   - For large reference genomes, ensure adequate memory is available

4. **Inconsistent Results After Adding Samples**
   - If results appear inconsistent after adding new samples, check that the same reference genome and parameters are being used
   - In rare cases, manually removing specific output files might be necessary to force reprocessing


## In developmnent

 - Working on wrapping entire pipeline into Singularity image
 - Intergrating snpEFF and snpSift scripts into main pipeline

## Citation

If you use this pipeline in your research, please cite:
[Insert citation information here]

## Contact
Please open an issue if you run into any problems