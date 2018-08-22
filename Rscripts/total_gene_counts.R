source('/home/ev250/Cincinatti/Functions/various.R')


## get args from snakemake
ebg <-  snakemake@input[[1]]

path <-  dirname(snakemake@input[[2]])
bam.name <- basename(snakemake@input[[2]])

mode <- snakemake@params[['mode']]
ignore.strand <- as.logical(snakemake@params[['ignore_strand']])

out <- snakemake@output[[1]]

## transform args to feed function

sample <- basename(path)
reads <- basename(dirname(path))

ebg <- readRDS(ebg)

singleEnd <- ifelse(reads=="Paired", FALSE, TRUE)

## Prepare matrix of counts per gene:

##cat(class(ebg), path, bam.name, mode, singleEnd, ignore.strand, out)

counts <- counts_sample_sub(ebg,path,bam.name, mode,singleEnd, ignore.strand)

names(counts)[2:ncol(counts)] <- sample

write.table(counts, file=out , row.names=F)
