# Work in progress...

# About
Workflow and code used in [Lewis R. et al 2025](https://www.biorxiv.org/content/10.1101/2024.12.23.628302v1.full) manuscript

# Table of contents
- [Abstract](#Abstract)
- [Reference genomes and annotations](#Reference-genomes-and-annotations)
- [Computational platform](#Computational-platform)
- [Requirements](#Requirements)
- [Indexing](#Indexing)
- [Read processing and alignment](#Read-processing-and-alignment)
    - [RNA-seq](#RNA-seq)

# Abstract
In most eukaryotic cells, euchromatin is localized in the nuclear interior, whereas heterochromatin is enriched at the nuclear envelope (NE). This conventional chromatin organization is established by heterochromatin tethering to the NE, however its importance for cellular homeostasis is largely unexplored. Peripheral heterochromatin localization relies on redundant NE-tethering systems. One tether is constituted by the lamin B receptor (LBR) in mammals, but the enigmatic nature of the other tethers has hampered functional analyses. Here we demonstrate that the downregulation of abundant, ubiquitous NE proteins can induce the global detachment of heterochromatin from the NE. Among these factors, we identify LBR and LAP2 as major players in bulk heterochromatin attachment to the NE in pluripotent and differentiated mammalian cells. Their loss leads to repositioning of heterochromatin to the nuclear interior, changes in chromatin accessibility, deregulation of gene expression including activation of antiviral innate immunity, and defects in cell fate determination.

# Reference genomes and annotations
The reference genomes for human ([GRCh38](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz); [annotation version 108](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz)) and mouse ([GRCm39](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz); [annotation version 112](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)) were obtained from ENSEMBL.

The following TEtranscripts annotation files were used to identify deregulated TE familes:
- [Humans](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm39_Ensembl_rmsk_TE.gtf.gz)
- [Mouse](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_Ensembl_rmsk_TE.gtf.gz)

# Computational platform
Genome indexing, read alignment and major parts of the processing for next generation sequencing data were done on [Euler HPC cluster](https://scicomp.ethz.ch/wiki/Euler) at ETH, Zurich.

The modules in the code refers to [applications and libraries](https://scicomp.ethz.ch/wiki/Euler_applications_and_libraries_ubuntu) available on Euler.

# Requirements

- STAR
- Bowtie2
- Samtools
- Cutadapt
- R

R libraies:
- Tidyverse
- DESeq2


# Indexing

Following is an example code showing indexing of the human genome using STAR. For mouse sample, please replace human genome and annotations with their mouse counterparts from above.


```bash

# Load modules
module load stack/2024-06
module load python/3.11.6
module load star/2.7.10b
module load samtools/1.17
module load pigz/2.7-oktqzxd

# Directories
dir_in="/path/to/dir"
gen_in="/path/to/dir"
THREADS=128

# Safety first
set -e

# Increase ulimit for STAR
ulimit -n 10240

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Unzipping Genome and Annotation"
echo "-------------------------"

# Unzip for STAR
pigz -p "$THREADS" -d "$gen_in"/*.gz

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Done Unzipping Genome and Annotation"
echo "-------------------------"

echo "-------------------------"
echo "Indexing Genome using STAR"
echo "-------------------------"

# Build Index
STAR \
--runThreadN "$THREADS" \
--runMode genomeGenerate \
--genomeDir "$gen_in"/ \
--genomeFastaFiles "$gen_in"/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--sjdbGTFfile "$gen_in"/Homo_sapiens.GRCh38.108.gtf \
--sjdbOverhang 99

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Done indexing"
echo "-------------------------"


```

# Read processing and alignment

## RNA-seq

Following is an example code showing processing of the HCT116 RNA-seq samples. For mESC samples, please replace human genome and annotations with their mouse counterparts from above.

The code does the following:
- Removes Illumina adaptors from the reads and keep the reads with `-q 25`
- Map the reads to the genome
- Index the coordinate sorted BAM file
- Count reads mapped to `exons` and `biotypes` and export them to their respective files.
- Performs QC on original fastq, processed fastq and BAM files
- Map the reads to genome with `--winAnchorMultimapNmax 200` and `--outFilterMultimapNmax 100` to generate BAM files for TEtranscripts analysis

Note: `featureCounts` crashes when `-T` > 64

```bash

# Load modules
module load stack/2024-06
module load python/3.11.6
module load star/2.7.10b
module load openjdk/17.0.8.1_1
module load fastqc/0.12.1
module load samtools/1.17
module load subread/2.0.6
module load py-dnaio/0.10.0-l2umcaf
module load py-xopen/1.6.0-o55adyx
module load py-cutadapt/4.4-jfcyzb5
module load pigz/2.7-oktqzxd

# Directories
dir_in="/path/to/dir"
gen_in="/path/to/dir"
THREADS=128

# Safety first
set -e

# Increase ulimit for STAR
ulimit -n 10240

# QC reads
for file in "$dir_in"/*_R1.fastq.gz; do
    base=$(basename "$file" _R1.fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cutadapt working on "${base}""
    echo "-------------------------"

    cutadapt \
    -j "$THREADS" \
    -q 25 \
    -m 25 \
    --poly-a \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o "$dir_in"/${base}_qc_R1.fastq.gz \
    -p "$dir_in"/${base}_qc_R2.fastq.gz \
    "$dir_in"/${base}_R1.fastq.gz \
    "$dir_in"/${base}_R2.fastq.gz \
    1> "$dir_in"/${base}_qc_cutadapt_log.txt

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cutadapt done processing $base"
    echo "-------------------------"

done


# Map with STAR
for file in "$dir_in"/*qc_R1.fastq.gz; do

    base=$(basename "$file" _qc_R1.fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "STAR aligning $base"
    echo "-------------------------"

    STAR \
        --runThreadN "$THREADS" \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir "$gen_in"/ \
        --readFilesCommand gunzip -c \
        --readFilesIn "$dir_in"/${base}_qc_R1.fastq.gz "$dir_in"/${base}_qc_R2.fastq.gz \
        --outFileNamePrefix "$dir_in"/${base}_

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "STAR done aligning $base"
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool indexing: "${base}_Aligned.sortedByCoord.out.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_Aligned.sortedByCoord.out.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE indexing: "${base}_Aligned.sortedByCoord.out.bam""
    echo "-------------------------"

done

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "featureCounts initiating"
echo "-------------------------"

# Exon Counts
featureCounts \
    -p \
    -T 64 \
    -t exon \
    -g gene_id \
    -a "$gen_in"/Homo_sapiens.GRCh38.108.gtf \
    -o "$dir_in"/HCT116_exon_counts_rawfile_tmp.txt \
    "$dir_in"/*Aligned.sortedByCoord.out.bam 2> "$dir_in"/Featurecount_exon_count_log.txt

# Biotypes
featureCounts \
    -p \
    -T 64 \
    -t exon \
    -g gene_biotype \
    -a "$gen_in"/Homo_sapiens.GRCh38.108.gtf \
    -o "$dir_in"/HCT116_biotype_counts_rawfile_tmp.txt \
    "$dir_in"/*Aligned.sortedByCoord.out.bam 2> "$dir_in"/Featurecount_biotype_counts_log.txt

# Remove first row
awk '(NR>1)' "$dir_in"/HCT116_exon_counts_rawfile_tmp.txt > "$dir_in"/HCT116_exon_counts_rawfile.txt
awk '(NR>1)' "$dir_in"/HCT116_biotype_counts_rawfile_tmp.txt > "$dir_in"/HCT116_biotype_counts_rawfile.txt

# For HCT
cut -f 1,7-41 "$dir_in"/HCT116_exon_counts_rawfile.txt > "$dir_in"/HCT116_exon_counts_processed.txt
cut -f 1,7-41 "$dir_in"/HCT116_biotype_counts_rawfile.txt > "$dir_in"/HCT116_biotype_counts_processed.txt

rm "$dir_in"/HCT116_exon_counts_rawfile_tmp.txt
rm "$dir_in"/HCT116_biotype_counts_rawfile_tmp.txt

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "featureCounts done!"
echo "-------------------------"


echo "-------------------------"
echo "FastQC working on FASTQ files"
echo "-------------------------"

# Make fastqc directory
mkdir "$dir_in"/fastqc

# FastQC on fastq files
fastqc \
-t "$THREADS" \
-o "$dir_in"/fastqc \
-f fastq "$dir_in"/*.fastq.gz

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "FastQC working on BAM files"
echo "-------------------------"

# FastQC on mapped files
fastqc \
-t "$THREADS" \
-o "$dir_in"/fastqc \
-f bam "$dir_in"/*.bam

echo "-------------------------"
echo "FastQC done!"
echo "-------------------------"


date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Starting mapping for TEtranscripts"
echo "-------------------------"

# Map with STAR for TEtranscripts
for file in "$dir_in"/*qc_R1.fastq.gz; do

    base=$(basename "$file" _qc_R1.fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "STAR aligning "${base}""
    echo "-------------------------"

    STAR \
        --runThreadN "$THREADS" \
        --outSAMtype BAM Unsorted \
        --sjdbOverhang 99 \
        --winAnchorMultimapNmax 200 \
        --outFilterMultimapNmax 100 \
        --genomeDir "$gen_in"/ \
        --readFilesCommand gunzip -c \
        --readFilesIn "$dir_in"/${base}_qc_R1.fastq.gz "$dir_in"/${base}_qc_R2.fastq.gz \
        --outFileNamePrefix "$dir_in"/${base}_forTE_

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "STAR DONE aligning "${base}""
    echo "-------------------------"

done

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Done!"
echo "-------------------------"

```



