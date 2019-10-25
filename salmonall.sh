#!/bin/bash
for f1 in *_R1.fastq.gz
do
        f2=${f1%%_R1.fastq.gz}"_R2.fastq.gz"
        /usr/home/b/485/tuj96866/salmon_0.99.0_beta2_linux_x86_64/bin/salmon quant --seqBias -l A -p 8 -i $1 -1 $f1 -2 $f2 --validateMappings -o "counts-$f1" 
done