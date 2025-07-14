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
    - [ATAC-seq](#ATAC-seq)
    - [ChIP-seq](#ChIP-seq)
- [Citation](#Citation)

# Abstract
In most eukaryotic cells, euchromatin is localized in the nuclear interior, whereas heterochromatin is enriched at the nuclear envelope (NE). This conventional chromatin organization is established by heterochromatin tethering to the NE, however its importance for cellular homeostasis is largely unexplored. Peripheral heterochromatin localization relies on redundant NE-tethering systems. One tether is constituted by the lamin B receptor (LBR) in mammals, but the enigmatic nature of the other tethers has hampered functional analyses. Here we demonstrate that the downregulation of abundant, ubiquitous NE proteins can induce the global detachment of heterochromatin from the NE. Among these factors, we identify LBR and LAP2 as major players in bulk heterochromatin attachment to the NE in pluripotent and differentiated mammalian cells. Their loss leads to repositioning of heterochromatin to the nuclear interior, changes in chromatin accessibility, deregulation of gene expression including activation of antiviral innate immunity, and defects in cell fate determination.

# Reference genomes and annotations
The reference genomes for human ([GRCh38](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz); [annotation version 108](https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz)) and mouse ([GRCm39](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz); [annotation version 112](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)) were obtained from ENSEMBL.

The following TEtranscripts annotation files were used to identify deregulated TE familes:
- [Humans](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCm39_Ensembl_rmsk_TE.gtf.gz)
- [Mouse](https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_Ensembl_rmsk_TE.gtf.gz)

# Computational platform
Genome indexing, read alignment and major parts of the processing for next generation sequencing data were done on [Euler HPC cluster](https://scicomp.ethz.ch/wiki/Euler) at ETH, Zurich.

The modules in the codes below refers to [applications and libraries](https://scicomp.ethz.ch/wiki/Euler_applications_and_libraries_ubuntu) available on Euler.

# Requirements

- STAR
- Bowtie2
- Samtools
- Cutadapt
- Subread
- Bedtools
- Sambamba
- MACS3
- R

R libraies:
- Tidyverse
- DESeq2


# Indexing

Following is an example code showing indexing of the human genome. For mouse sample, please replace human genome and annotations with their mouse counterparts from above.


```bash

# Load modules
module load stack/2024-06
module load python/3.11.6
module load star/2.7.10b
module load bowtie2/2.5.1-u2j3omo
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

# Inde genome for bowtie2
date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Indexing Genome using Bowtie2"
echo "-------------------------"

bowtie2-build \
--threads "$THREADS" \
"$gen_in"/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
"$gen_in"/Homo_sapiens.GRCh38.dna.primary_assembly

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Indexing done!"
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

cut -f 1,7- "$dir_in"/HCT116_exon_counts_rawfile.txt > "$dir_in"/HCT116_exon_counts_processed.txt
cut -f 1,7- "$dir_in"/HCT116_biotype_counts_rawfile.txt > "$dir_in"/HCT116_biotype_counts_processed.txt

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

## ATAC-seq

The code does the following:
- Removes Illumina adaptors from the reads and keep the reads with `-q 25`
- Map the reads to the genome, sort by coordinates and index them
- Remove mitochondrial reads
- Keep only uniquely mapped reads
- Checks fragment length distribution (QC)
- Select reads with 100 nt length and shift the 5' and 3' ends to account for transposon dimerization before insertion
- Call peaks on above file using MACS3
- Filter out black listed regions from the above file.

```bash

module load stack/2024-06
module load python/3.11.6
module load bowtie2/2.5.1-u2j3omo
module load samtools/1.17
module load subread/2.0.6
module load py-dnaio/0.10.0-l2umcaf
module load py-xopen/1.6.0-o55adyx
module load py-cutadapt/4.4-jfcyzb5
module load openjdk/17.0.8.1_1
module load fastqc/0.12.1
module load r/4.4.0
module load xz/5.4.1-hrbyxav
module load bzip2/1.0.8-kxrmhba
module load curl/8.4.0-s6dtj75
module load libxml2/2.10.3-xbqziof
module load libiconv/1.17-uiaqkl2
module load pigz/2.7-oktqzxd
module load bedtools2/2.31.0


# Directories
dir_in="/path/to/dir"
gen_in="/path/to/dir"
smabamba_dir="/path/to/dir/sambamba_v1.0.1"
deeptools_dir="/path/to/dir/"
picard_dir="/path/to/dir/picard/build/libs"
tools_dir="/path/to/dir"
THREADS=128

# Safety first
set -e


# Trim reads
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
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -o "$dir_in"/${base}_qc_R1.fastq.gz \
    -p "$dir_in"/${base}_qc_R2.fastq.gz \
    "$dir_in"/${base}_R1.fastq.gz \
    "$dir_in"/${base}_R2.fastq.gz \
    1> "$dir_in"/${base}_qc_cutadapt_log.txt

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cutadapt done processing "${base}""
    echo "-------------------------"

done


# Map with bowtie2
for file in "$dir_in"/*qc_R1.fastq.gz; do

    base=$(basename "$file" _qc_R1.fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bowtie2 aligning "${base}" for MACS3"
    echo "-------------------------"

    (
        bowtie2 \
        -t \
        -q \
        -p 64 \
        --very-sensitive \
        -x "$gen_in"/Homo_sapiens.GRCh38.dna.primary_assembly \
        -1 "$dir_in"/${base}_qc_R1.fastq.gz \
        -2 "$dir_in"/${base}_qc_R2.fastq.gz | \
        samtools view \
        -@ 64 -h -b -S - > "$dir_in"/${base}_macs3_genome_unsorted.bam
    ) 2> "$dir_in"/${base}_macs3_genome_map_log.txt

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bowtie2 DONE aligning "${base}" for MACS3"
    echo "-------------------------"

    

    echo "-------------------------"
    echo "Samtool sorting: "${base}_macs3_genome_unsorted.bam""
    echo "-------------------------"

    # Sort the BAM file
    samtools sort \
    -@ "$THREADS" \
    -o "$dir_in"/${base}_macs3_genome_sorted.bam "$dir_in"/${base}_macs3_genome_unsorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE sorting: "${base}_macs3_genome_unsorted.bam""
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_macs3_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE indexing: "${base}_macs3_genome_sorted.bam""
    echo "-------------------------"

    echo "-------------------------"
    echo "Cleaning up: unsorted bams"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_genome_unsorted.bam

done

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Bowtie2 and SAMtools done!"
echo "-------------------------"

echo "-------------------------"
echo "Processing BAMs for MACS3!"
echo "-------------------------"


# Remove MT reads
for file in "$dir_in"/*_macs3_genome_sorted.bam; do

    base=$(basename "$file" _macs3_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "SAMtools working on "${base}""
    echo "-------------------------"

    samtools idxstats \
    "$dir_in"/${base}_macs3_genome_sorted.bam | \
    cut -f1 | \
    grep -v MT | \
    xargs samtools view \
    --threads 54 \
    -b \
    "$dir_in"/${base}_macs3_genome_sorted.bam > "$dir_in"/${base}_macs3_nomt_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "SAMtools done processing $base"
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_nomt_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ 64 \
    "$dir_in"/${base}_macs3_nomt_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cleaning up: bams with MT reads"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_genome_sorted.bam*

done

# Remove duplicates and keep only uniquely mapped reads
for file in "$dir_in"/*_macs3_nomt_genome_sorted.bam; do

    base=$(basename "$file" _macs3_nomt_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Sambamba working on "${base}""
    echo "-------------------------"

    "$smabamba_dir"/sambamba_v1.0.1 view \
    -h \
    -t "$THREADS" \
    -f bam \
    -F "[XS] == null and not unmapped and not duplicate" \
    "$dir_in"/${base}_macs3_nomt_genome_sorted.bam > "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam

    echo "-------------------------"
    echo "Sambamba done processing "${base}""
    echo "-------------------------"

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_nomt_unique_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cleaning up: bams with duplicate reads"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_nomt_genome_sorted.bam*

done

# Keep only nucleosome free regions (< 100 bp reads) and shift reads

source "$deeptools_dir"/deeptools/bin/activate

for file in "$dir_in"/*_macs3_nomt_unique_genome_sorted.bam; do

    base=$(basename "$file" _macs3_nomt_unique_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Deeptools working on "${base}""
    echo "-------------------------"

    alignmentSieve \
    --numberOfProcessors "$THREADS" \
    --minMappingQuality 25 \
    --maxFragmentLength 100 \
    --ATACshift \
    --BED \
    --filterMetrics log.txt \
    --bam "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam \
    --outFile "$dir_in"/${base}_macs3_input.bed  

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Deeptools done processing "${base}""
    echo "-------------------------"


done

deactivate


# Fragment length distribution
for file in "$dir_in"/*_macs3_nomt_unique_genome_sorted.bam; do

    base=$(basename "$file" _macs3_nomt_unique_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Picard working on "${base}""
    echo "-------------------------"

    java \
    -Xmx64G \
    -jar "$picard_dir"/picard.jar \
    CollectInsertSizeMetrics \
    -M 0.5 \
    -I "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam \
    -O "$dir_in"/${base}_frag_len_stats.txt \
    -H "$dir_in"/${base}_frag_len.pdf

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Picard done processing "${base}""
    echo "-------------------------"


done

# Call peaks

source "$tools_dir"/macs3/bin/activate

for file in "$dir_in"/*_macs3_input.bed; do

    base=$(basename "$file" _macs3_input.bed)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "MACS3 working on "${base}""
    echo "-------------------------"

    macs3 \
    callpeak \
    --call-summits \
    --nomodel \
    --nolambda \
    --gsize hs \
    --keep-dup all \
    --format BEDPE \
    --qvalue 0.05 \
    --outdir "$dir_in"/ \
    --treatment "$file" \
    --name "$dir_in"/${base}_macs3 \
    2> "$dir_in"/${base}_macs3_log.txt

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "MACS3 done processing "${base}""
    echo "-------------------------"


done

deactivate

# Filter-out blaclisted regions

for file in "$dir_in"/*_all_macs3_peaks.narrowPeak; do

    base=$(basename "$file" _all_macs3_peaks.narrowPeak)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bedtools working on "${base}""
    echo "-------------------------"

    bedtools intersect \
    -v \
    -a "$file" \
    -b "$gen_in"/hg38_blacklist_v2_ps.bed \
    > "$dir_in"/${base}_all_macs3_filtered_peaks.bed 

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bedtools done processing "${base}""
    echo "-------------------------"

done

```

## ChIP-seq

The code does the following:
- Removes Illumina adaptors from the reads and keep the reads with `-q 25`
- Map the reads to the genome, sort by coordinates and index them
- Remove mitochondrial reads
- Keep only uniquely mapped reads
- Remove blacklisted region from the BAM files
- Sort the above BAM files based on coordinates and index them


```bash

module load stack/2024-06
module load python/3.11.6
module load bowtie2/2.5.1-u2j3omo
module load samtools/1.17
module load subread/2.0.6
module load py-dnaio/0.10.0-l2umcaf
module load py-xopen/1.6.0-o55adyx
module load py-cutadapt/4.4-jfcyzb5
module load openjdk/17.0.8.1_1
module load fastqc/0.12.1
module load r/4.4.0
module load xz/5.4.1-hrbyxav
module load bzip2/1.0.8-kxrmhba
module load curl/8.4.0-s6dtj75
module load libxml2/2.10.3-xbqziof
module load libiconv/1.17-uiaqkl2
module load bedtools2/2.31.0
module load pigz/2.7-oktqzxd

# Directories
dir_in="/path/to/dir"
gen_in="/path/to/dir"
smabamba_dir="/path/to/dir/sambamba_v1.0.1"
tools_dir="/path/to/dir/"
THREADS=128

# Safety first
set -e


# Trim reads
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
    echo "Cutadapt done processing "${base}""
    echo "-------------------------"

done


# Map with bowtie2
for file in "$dir_in"/*qc_R1.fastq.gz; do

    base=$(basename "$file" _qc_R1.fastq.gz)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bowtie2 aligning "${base}" for MACS3"
    echo "-------------------------"

    (
        bowtie2 \
        -t \
        -q \
        -p 64 \
        --very-sensitive \
        -x "$gen_in"/Homo_sapiens.GRCh38.dna.primary_assembly \
        -1 "$dir_in"/${base}_qc_R1.fastq.gz \
        -2 "$dir_in"/${base}_qc_R2.fastq.gz | \
        samtools view \
        -@ 64 -h -b -S - > "$dir_in"/${base}_macs3_genome_unsorted.bam
    ) 2> "$dir_in"/${base}_macs3_genome_map_log.txt

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bowtie2 DONE aligning "${base}" for MACS3"
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool sorting: "${base}_macs3_genome_unsorted.bam""
    echo "-------------------------"

    # Sort the BAM file
    samtools sort \
    -@ "$THREADS" \
    -o "$dir_in"/${base}_macs3_genome_sorted.bam "$dir_in"/${base}_macs3_genome_unsorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE sorting: "${base}_macs3_genome_unsorted.bam""
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_macs3_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE indexing: "${base}_macs3_genome_sorted.bam""
    echo "-------------------------"

    echo "-------------------------"
    echo "Cleaning up: unsorted bams"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_genome_unsorted.bam

done

date +"%d-%m-%Y %T"
echo "-------------------------"
echo "Bowtie2 and SAMtools done!"
echo "-------------------------"

echo "-------------------------"
echo "Processing BAMs for MACS3!"
echo "-------------------------"


# Remove MT reads
for file in "$dir_in"/*_macs3_genome_sorted.bam; do

    base=$(basename "$file" _macs3_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "SAMtools working on "${base}""
    echo "-------------------------"

    samtools idxstats \
    "$dir_in"/${base}_macs3_genome_sorted.bam | \
    cut -f1 | \
    grep -v MT | \
    xargs samtools view \
    --threads "$THREADS" \
    -b \
    "$dir_in"/${base}_macs3_genome_sorted.bam > "$dir_in"/${base}_macs3_nomt_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "SAMtools done processing $base"
    echo "-------------------------"

    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_nomt_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_macs3_nomt_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cleaning up: bams with MT reads"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_genome_sorted.bam*

done

# Remove duplicates and keep only uniquely mapped reads
for file in "$dir_in"/*_macs3_nomt_genome_sorted.bam; do

    base=$(basename "$file" _macs3_nomt_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Sambamba working on "${base}""
    echo "-------------------------"

    "$smabamba_dir"/sambamba_v1.0.1 view \
    -h \
    -t "$THREADS" \
    -f bam \
    -F "[XS] == null and not unmapped and not duplicate" \
    "$dir_in"/${base}_macs3_nomt_genome_sorted.bam > "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam

    echo "-------------------------"
    echo "Sambamba done processing "${base}""
    echo "-------------------------"

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool indexing: "${base}_macs3_nomt_unique_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Cleaning up: bams with duplicate reads"
    echo "-------------------------"

    rm "$dir_in"/${base}_macs3_nomt_genome_sorted.bam*

done

# Remove blacklisted regions
for file in "$dir_in"/*_macs3_nomt_unique_genome_sorted.bam; do

    base=$(basename "$file" _macs3_nomt_unique_genome_sorted.bam)

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bedtools working on "${base}""
    echo "-------------------------"

    bedtools intersect \
    -v \
    -abam "$dir_in"/${base}_macs3_nomt_unique_genome_sorted.bam \
    -b "$gen_in"/hg38_blacklist_v2_ps.bed | \
    samtools sort \
    -@ "$THREADS" \
    -o "$dir_in"/${base}_nomt_noblklstd_unique_genome_sorted.bam -

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Bedtools DONE working on "${base}""
    echo "-------------------------"

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool indexing: "${base}_nomt_noblklstd_unique_genome_sorted.bam""
    echo "-------------------------"

    # Index the BAM file
    samtools index \
    -@ "$THREADS" \
    "$dir_in"/${base}_nomt_noblklstd_unique_genome_sorted.bam

    date +"%d-%m-%Y %T"
    echo "-------------------------"
    echo "Samtool DONE indexing: "${base}_nomt_noblklstd_unique_genome_sorted.bam""
    echo "-------------------------"


done


```


# Citation
LBR and LAP2 mediate heterochromatin tethering to the nuclear periphery to preserve genome homeostasis
*Renard Lewis, Virginia Sinigiani, Krisztian Koos, Cristiana Bersaglieri, Caroline Ashiono, Raffaella Santoro, Constance Ciaudo, Peter Horvath, Puneet Sharma, Ulrike Kutay*
bioRxiv 2024.12.23.628302; doi: [https://doi.org/10.1101/2024.12.23.628302](https://doi.org/10.1101/2024.12.23.628302)