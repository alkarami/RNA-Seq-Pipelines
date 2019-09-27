## Welcome to my tutorial on DESeq2! A lot of this is derived from the official documentation, but
## some functions in the original paper have been deprecated - so I thought, why not write a tutorial
## that incorporates the working functions and show you how it would be applied to an actual data set?

## I used these exact lines of code to interrogate difrerential gene expression in the cell lines 
## my lab works on. 

## Your first step is to define your genes and gene loci. This is done by obtaining and uploading a 
## .gtf file or "gene transfer file"
## Get a .gtf for your organism of interest. As I was working on human cell lines, I got the latest 
## human .gtf from GENCODE. 

library( "GenomicFeatures" )
hse <- makeTxDbFromGFF( "gencode.v31.annotation.gtf", format="gtf" )
exonsByGene <- exonsBy( hse, by="gene" )

## Now you should define the directory in which you will store all your .bam files. 
## In this experiment, my files were all in a folder called "Reshu_bams" in my directory.

fls <- list.files( "Reshu_bams", pattern="bam$", full=TRUE )

## Only process 100000 reads at once. You can change this to lower or higher given your processing
## power. With my 16G system, 100000 reads was fine. Could probably have gone higher, but whatever.

library( "Rsamtools" )
bamLst <- BamFileList( fls, yieldSize=100000 )

## Summarize overlaps. summarizeOverlaps is a core function for what we're trying to do as it 
## summarizes the overlaps between the bams' read mapping and your defined list of genes. This is the 
## gas guzzler of your code. 

library( "GenomicAlignments" )
se <- summarizeOverlaps( exonsByGene, bamLst,
                         mode="Union",
                         singleEnd=TRUE,
                         ignore.strand=FALSE,
                         fragments=FALSE )

## Finally, now you can load up DESeq2 and make a DESeq2 dataset with the summarized experiment.
## Upload your metadata so that you can easily define what your samples are and what groups they 
## all belong to. The .csv that I had had two columns: samples and group. Samples had the sample names
## and group had the timepoint from which the samples were derived. You can define these columns
## however you like, just make sure your "samples" column match your .bam filenames.
## The "group" defines the technical replicates. The metadata file is useful because now you can 
## insert any information you want pertaining to each sample and compare and contrast accordingly

library( "DESeq2" )
sampleInfo <- read.csv( "reshu_meta.csv" )
seIdx <- match(colnames(se), sampleInfo$samples)
colData(se) <- cbind(colData(se),sampleInfo[ seIdx, ])
ddsFull <- DESeqDataSet( se, design = ~group)

## Run DESeq. Shouldn't take too long. 
dds <- DESeq(ddsFull)

## Congrats! Now you can extract the results. Follow the lines below to get the normalized counts 
## for each gene for each sample

dds <- estimateSizeFactors(dds)
counts(dds, normalized=TRUE)

## You can also get results for fold changes, which was what I was most interested in. 
## The first argument is the column from which you're getting the contrast from, and 
## the 2nd and 3rd arguments define the numerator and denominator, respectively. 

res3v1 <- results( dds, contrast = c("group", "3d", "cont") )
res7v3 <- results( dds, contrast = c("group", "7d", "3d") )
res7v1 <- results( dds, contrast = c("group", "7d", "cont") )


## Annotate with gene names! Check out your results - the rownames reveal each and every loci 
## but for downstream purposes you will want to get actual gene names. Match rownames with gene names
## by accessing a database like org.Hs.eg.db

## Some gene builds (like the one I used in the example) will use ENSEMBL gene name with the accession
## numbers (I.e. ENSG90022123.1 instead of ENSG90022123). That's not good for downstream gene search
## efforts, as the '.#' part of the string ruins the search. Here I use a package to take care of 
## that (because I'm bad at R). 

library(stringr)
gene_ids <- str_replace(row.names(dds),
                        pattern = ".[0-9]+$",
                        replacement = "")
library("AnnotationDbi")
library("org.Hs.eg.db")
symbols <- mapIds(org.Hs.eg.db,
                     keys=gene_ids,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

## Cool, now you can assign gene names (gene symbols) as an extra column in the data frame. 

res3v1$symbol <- symbols
res7v3$symbol <- symbols
res7v1$symbol <- symbols

## Print results into a csv if you want!
write.csv( as.data.frame(res3v1), file="results3v1_DESeq2.csv" )
write.csv( as.data.frame(res7v3), file="results7v3_DESeq2.csv" )
write.csv( as.data.frame(res7v1), file="results7v1_DESeq2.csv" )
