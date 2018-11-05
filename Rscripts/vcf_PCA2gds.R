library(gdsfmt)
library(SNPRelate)
library(parallel)

## Reads data from vcf files and convert into SNP GDS Format, calls snpgdsVCF2GDS from SNPRelate package

len=length(snakemake@input)

vcf.in = list(snakemake@input[1:(len/2)], snakemake@input[(len/2+1):len])
## convert each vcf into gds and then merge them: merge gds didnt work

##peac.vcf <- snakemake@input[1:(len/2)]
##rp.vcf <- snakemake@input[(len/2+1):len]

mclapply(1:2 ,function(i)  snpgdsVCF2GDS(vcf.fn=unlist(vcf.in[[i]]), method=snakemake@params[['method']],out.fn=snakemake@output[[i]]), mc.cores = snakemake@threads)


## ## peac
## mclapply( 1:length(peac.vcf) ,
##          function(i)  snpgdsVCF2GDS(vcf.fn=peac.vcf[[i]], method=snakemake@params[['method']], out.fn=paste0(snakemake@params[['prefix']],"PEAC_", i,".gds" )),
##          mc.cores = snakemake@threads)

## ## rp
## mclapply( 1:length(rp.vcf) ,function(i)  snpgdsVCF2GDS(vcf.fn=rp.vcf[[i]], method=snakemake@params[['method']],out.fn=paste0(snakemake@params[['prefix']],"RP_", i,".gds" )), mc.cores = snakemake@threads)

## ## merge

## peac <- list.files(path=snakemake@params[['prefix']], pattern="PEAC_[0-9]+.gds", full.names=TRUE)

## rp <- list.files(path=snakemake@params[['prefix']], pattern="RP_[0-9]+.gds", full.names=TRUE)

## snpgdsCombineGeno(gds.fn=peac, out.fn=snakemake@output[['peac']])
## snpgdsCombineGeno(gds.fn=rp, out.fn=snakemake@output[['rp']])
 
## ## remove temp files

## file.remove(peac)
## file.remove(rp)

