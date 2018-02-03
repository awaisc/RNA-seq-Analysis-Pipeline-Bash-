#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=10:00:00
#SBATCH --mem=32GB
#SBATCH -o /fast/users/a1649239/GeneAnalysis/AddrenalGlandNasa/StevesScriptOutput/runPipeline_%j.out
#SBATCH -e /fast/users/a1649239/GeneAnalysis/AddrenalGlandNasa/StevesScriptOutput/runPipeline_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## Full run of the Hypoxia data.

## Modules
module load fastqc/0.11.4
module load HISAT2/2.1.0-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load AdapterRemoval/2.2.1-foss-2016b
module load GCC/5.4.0-2.26

## Biohub/local
sambamba=/data/biohub/local/sambamba_v0.6.6
featureCounts=/data/biohub/local/subread-1.5.2-Linux-x86_64/bin/featureCounts

## Directories
ROOT=/fast/users/a1649239/GeneAnalysis/TestNasaPairedEnd
REFS=/data/biohub/Refs/Mus_musculus/ensembl85
MERGEDATA=${ROOT}/01_mergedData/fastq
TRIMDATA=${ROOT}/02_trimmedData
ALIGNDATA=${ROOT}/03_alignedData

## Making directories for Trimmed data
mkdir -p ${TRIMDATA}
mkdir -p ${TRIMDATA}/fastq
mkdir -p ${TRIMDATA}/fastqMerged
mkdir -p ${TRIMDATA}/FastQC

## Making directories for zebrafish genome alignment
mkdir -p ${ALIGNDATA}
mkdir -p ${ALIGNDATA}/sams
mkdir -p ${ALIGNDATA}/logs
mkdir -p ${ALIGNDATA}/bams
mkdir -p ${ALIGNDATA}/duplicatesRemoved
mkdir -p ${ALIGNDATA}/FastQC
mkdir -p ${ALIGNDATA}/featureCounts


## Data files
z10_hisat2index=${REFS}/HiSat2/grcm38/genome
z10_gtf=${REFS}/gencode.vM11.primary_assembly.annotation.gtf

## Cores
CORES=32


module load picard/2.2.4-Java-1.8.0_71

java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar


for BAMFILES in /fast/users/a1649239/GeneAnalysis/TestNasaPairedEnd/03_alignedData/bams/*.bam
do

java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar MarkDuplicates INPUT= ${BAMFILES} OUTPUT=${ALIGNDATA}/duplicatesRemoved/$(basename ${BAMFILES} trim1.fastq.gz_sorted.bam)Dedup_trim1.fastq.gz_sorted.bam METRICS_FILE= ${ALIGNDATA}/duplicatesRemoved/$(basename ${BAMFILES} trim1.fastq.gz_sorted.bam)Dedup_trim1.fastq.gz_sorted.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

done
