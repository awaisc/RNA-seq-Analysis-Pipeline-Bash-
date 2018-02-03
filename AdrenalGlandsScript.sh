#!/bin/bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH -o /fast/users/a1649239/GeneAnalysis/AddrenalGlandNasa/M23/runPipeline_%j.out
#SBATCH -e /fast/users/a1649239/GeneAnalysis/AddrenalGlandNasa/M23/runPipeline_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

## Full run of the Hypoxia data.

## Modules
module load fastqc/0.11.4
module load HISAT2/2.1.0-foss-2016b
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25
module load AdapterRemoval/2.2.1-foss-2016b
module load GCC/5.4.0-2.26
module load picard/2.2.4-Java-1.8.0_71

## Biohub/local
sambamba=/data/biohub/local/sambamba_v0.6.6
featureCounts=/data/biohub/local/subread-1.5.2-Linux-x86_64/bin/featureCounts

## Directories
ROOT=/fast/users/a1649239/GeneAnalysis/AddrenalGlandNasa/M23
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
z10_gtf=${REFS}/Mus_musculus.GRCm38.85.gtf


## Cores
CORES=32

##--------------------------------------------------------------------------------------------##
## Trimming the Merged data
##--------------------------------------------------------------------------------------------##

for FQGZ in ${ROOT}/*R1_001*.fastq.gz
do
echo "Trimming ${FQGZ}"

        AdapterRemoval \
        --file1 ${FQGZ} \
	--file2 ${FQGZ/R1_001/R2_001} \
        --output1 ${TRIMDATA}/fastq/$(basename $FQGZ _R1_001.fastq.gz)_trim1.fastq.gz \
	--output2 ${TRIMDATA}/fastq/$(basename $FQGZ _R1_001.fastq.gz)_trim2.fastq.gz \
        --gzip \
        --trimqualities --minquality 20 --minlength 35 \
        --threads ${CORES} \
        --qualitymax 50
done

##--------------------------------------------------------------------------------------------##
## Merging The Trimmed Data
##--------------------------------------------------------------------------------------------##


NumberOfUniqueBioSamples=$(ls ${TRIMDATA}/fastq/*trim1*fastq.gz | sed 's!.*/!!'  |   sed 's/-.*/ /' | uniq )

echo " Unique Biological samples dected ${NumberOfUniqueBioSamples}"

#Making a loop to isolate each sample name and number and merging these files together!

for uniqueSamples in $NumberOfUniqueBioSamples
do
echo "Merging lanes for bioloigcal Sample: $uniqueSamples"

cat ${TRIMDATA}/fastq/${uniqueSamples}*trim1* > ${TRIMDATA}/fastqMerged/${uniqueSamples}_R1_trim1.fastq.gz
cat ${TRIMDATA}/fastq/${uniqueSamples}*trim2* > ${TRIMDATA}/fastqMerged/${uniqueSamples}_R1_trim2.fastq.gz


echo "Done Merging files to ${TRIMDATA}/fastqMerged/${uniqueSamples}_R1_trim1.fastq.gz"
echo "Done Merging files to ${TRIMDATA}/fastqMerged/${uniqueSamples}_R1_trim2.fastq.gz"
done


##--------------------------------------------------------------------------------------------##
## FastQC on the merged Data
##--------------------------------------------------------------------------------------------##




fastqc -t ${CORES} -o ${TRIMDATA}/FastQC --noextract ${TRIMDATA}/fastqMerged/*


##--------------------------------------------------------------------------------------------##
## Aligning trimmed data to the zebrafish genome
##--------------------------------------------------------------------------------------------##

## Aligning, filtering and sorting
for FQGZ in ${TRIMDATA}/fastqMerged/*trim1.fastq.gz
do
    ## Alignment to ZF reference
    hisat2 -x ${z10_hisat2index} \
      -p ${CORES} \
      -1 ${FQGZ} \
      -2 ${FQGZ/trim1/trim2} \
      -S ${ALIGNDATA}/sams/$(basename ${FQGZ} _R1_001.fastq.gz).sam \
      &> ${ALIGNDATA}/logs/$(basename ${FQGZ} _R1_001.fastq.gz).log


    ## Convert SAM to BAM
    samtools view -bhS -q30 \
      ${ALIGNDATA}/sams/$(basename ${FQGZ} _R1_001.fastq.gz).sam \
      1> ${ALIGNDATA}/bams/$(basename ${FQGZ} _R1_001.fastq.gz).bam


    ## Sort & index bam files
    ${sambamba} sort -t ${CORES} \
      -o ${ALIGNDATA}/bams/$(basename ${FQGZ} _R1_001.fastq.gz)_sorted.bam \
      ${ALIGNDATA}/bams/$(basename ${FQGZ} _R1_001.fastq.gz).bam

    ## Remove the unsorted file
    rm ${ALIGNDATA}/bams/$(basename ${FQGZ} _R1_001.fastq.gz).bam

  done

# Fastqc
for BAM in ${ALIGNDATA}/bams/*_sorted.bam
  do
    fastqc -t ${CORES} -f bam_mapped -o ${ALIGNDATA}/FastQC --noextract ${BAM}
  done

  ##--------------------------------------------------------------------------------------------##
  ## Remove Duplicates Picard
  ##--------------------------------------------------------------------------------------------##

  for BAMFILES in /fast/users/a1649239/GeneAnalysis/TestNasaPairedEnd/03_alignedData/bams/*.bam
  do

  java -jar /apps/software/picard/2.2.4-Java-1.8.0_71/picard.jar MarkDuplicates INPUT= ${BAMFILES} OUTPUT=${ALIGNDATA}/duplicatesRemoved/$(basename ${BAMFILES} trim1.fastq.gz_sorted.bam)Dedup_trim1.fastq.gz_sorted.bam METRICS_FILE= ${ALIGNDATA}/duplicatesRemoved/$(basename ${BAMFILES} trim1.fastq.gz_sorted.bam)Dedup_trim1.fastq.gz_sorted.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT

  done


##--------------------------------------------------------------------------------------------##
## featureCounts
##--------------------------------------------------------------------------------------------##

## Feature Counts - obtaining all sorted bam files
sampleList=`find ${ALIGNDATA}/bams -name "*_sorted.bam" | tr '\n' ' '`

## Running featureCounts on the sorted bam files
${featureCounts} -Q 10 \
  -p \
  -T ${CORES} \
  -a ${z10_gtf} \
  -o ${ALIGNDATA}/featureCounts/counts.out ${sampleList}

## Storing the output in a single file
cut -f1,7- ${ALIGNDATA}/featureCounts/counts.out | \
sed 1d > ${ALIGNDATA}/featureCounts/genes.out
