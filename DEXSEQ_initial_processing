#!/bin/sh

##########################################################################
## A script template for submitting batch jobs. To submit a batch job, 
## please type
##
##    qsub myprog.sh
##
## Please note that anything after the first two characters "#$" on a line
## will be treated as a SUN Grid Engine command.
##########################################################################

## The following to run programs in the current working directory

#$ -cwd


## Specify a queue

#$ -q batchq


## The following two lines will send an email notification when the job is 
## Ended/Aborted/Suspended - Please replace "UserName" with your own username.

#$ -M aangeles
#$ -m eas

## Takes in argument from command line
file=$1

## Load modules
echo "$file is going to be processed"
module load trim_galore/0.3.1
module load fastqc/0.10.1
module load hisat/2.1.0
module load python3-cbrg

## Creates variable of file without _1.fastq and directory path
base=$(basename ${file/_1.fastq/});

## Trim adapters from raw sequencing data FASTQ files using Trim Galore
trim_galore --paired --retain_unpaired -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ${file} ${file/_1/_2} -o /t1-data/user/aangeles/fastaq_files/Test

## Wait until trimmed adapter files are generated before proceeding to next steps
while [ ! -e ${base}_2_val_2.fq ]; do sleep 10; done
while [ ! -e ${base}_1_val_1.fq ]; do sleep 10; done

## Generate FastQC reports from trimmed adapter files
fastqc -f fastq ${base}_1_val_1.fq ${base}_2_val_2.fq

## Align trimmed adapter files to genome using HISAT2 and generate BAM files
hisat2 -x hg38 -p 4 -1 ${base}_1_val_1.fq -2 ${base}_2_val_2.fq | samtools view -bSh > ${base}.bam

## Delete trimmed adapter files after HISAT2
rm ${base}_1_val_1.fq ${base}_2_val_2.fq
rm ${base}_1_unpaired_1.fq ${base}_2_unpaired_2.fq

## Generate count files from BAM files
python /t1-data/user/aangeles/fastaq_files/dexseq_count_mod.py /t1-data/user/aangeles/fastaq_files/Test/gencode_hg38_CACNA1C.gff ${base}.bam ${base}_CACNA1C.counts.txt; 
