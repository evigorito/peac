source('/home/ev250/Cincinatti/Functions/various.R')

## prepare exon by gene input file for calculating counts per gene 

ebg <- gtf2ebg(snakemake@input[[1]])

saveRDS(ebg, snakemake@output[[1]])




