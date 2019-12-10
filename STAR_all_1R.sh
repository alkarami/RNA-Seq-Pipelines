#!/bin/bash

## Arg 1: Name of the genome index; Arg 2: directory/name of the fasta file to make an index from 

## Make index directory
mkdir $1

## Make the index
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir gindex --genomeFastaFiles $2

## STAR every fastq.gz file
for f1 in *_R1_001.fastq.gz
do
	f3=${f1%%_R1_001.fastq.gz}
        nohup STAR --runThreadN 16 --genomeDir gindex --readFilesIn $f1 --readFilesCommand zcat --outFileNamePrefix "$f3" --outFilterMultimapNmax 1 --outSAMtype BAM SortedByCoordinate
done
