#!/bin/bash

## SALMON VERSION. END-TO-END FASTQ TO SALMON QUANT PROCESSOR. RUN THIS ON THE DIRECTORY WITH YOUR FOLDERS OF FASTQS
## GET CATLANES, FASTPALL, AND SALMONALL

## Unravel all subdirectories:
find . -type f -exec mv --backup=numbered {} . \; && find . -maxdepth 1 -type d -exec rm -r {} +

## Combine all fastq's from different lanes 
bash catlanes.sh 

## Run fastp on all fastq's (MAY NEED TO TINKER WITH THE ENDING OF FILENAME - IE WHAT COMES BEFORE FASTQ.GZ)
bash fastpall.sh

## Ready for salmon! Run salmon on all fastq's, specifying salmon index as the argument
nohup bash /usr/home/b/485/tuj96866/salmonall.sh /usr/home/b/485/tuj96866/mousein
