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
module load rna-star/2.6.1d
module load python3-cbrg

## Creates variable of file without _1.fastq and directory path
base=$(basename ${file/_1.fastq/});

## Trim adapters from raw sequencing data FASTQ files using Trim Galore
trim_galore --paired --retain_unpaired -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ${file} ${file/_1/_2} -o /t1-data/user/aangeles/fastaq_files/DESEQ

## Wait until trimmed adapter files are generated before proceeding to next steps
while [ ! -e ${base}_2_val_2.fq ]; do sleep 10; done
while [ ! -e ${base}_1_val_1.fq ]; do sleep 10; done

## Generate FastQC reports from trimmed adapter files
fastqc -f fastq ${base}_1_val_1.fq ${base}_2_val_2.fq

## Gunzip trimmed adapter files
gzip ${base}_1_val_1.fq
gzip ${base}_2_val_2.fq

## Align trimmed adapter files to genome using STAR
STAR --runThreadN 8 --genomeDir /t1-data/user/aangeles/fastaq_files/DESEQ --readFilesIn ${base}_1_val_1.fq.gz ${base}_2_val_2.fq.gz --outFileNamePrefix /t1-data/user/aangeles/fastaq_files/DESEQ/${base}_ --twopassMode Basic --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted

## Generate count files using HTSEQ
htseq-count --mode=intersection-strict --stranded=no --type=exon --format=bam --idattr=gene_id --quiet /t1-data/user/aangeles/fastaq_files/DESEQ/${base}_Aligned.out.bam /t1-data/user/aangeles/fastaq_files/new_annotation.gtf >${base}.txt

## Delete trimmed adapter files after STAR
rm ${base}_1_val_1.fq ${base}_2_val_2.fq
rm ${base}_1_unpaired_1.fq ${base}_2_unpaired_2.fq
