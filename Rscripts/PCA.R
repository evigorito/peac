library(gdsfmt)
library(SNPRelate)


#' Computes PCA from reference panel after pruning variants

## LD pruning
RPgenofile <- snpgdsOpen(snakemake@input[["RP"]])
snpset <- snpgdsLDpruning(RPgenofile, ld.threshold = as.numeric(snakemake@params[["ld"]]),
                          maf = as.numeric(snakemake@params[["maf"]]),
                          method = snakemake@params[["method"]],
                          verbose= TRUE)
##Pcs
pc <- snpgdsPCA(RPgenofile, snp.id=unlist(snpset), num.thread = snakemake@threads)
saveRDS(pc, file = snakemake@output[["PC"]])


## Matrix of loadings (rotation matrix)
snpL <-  snpgdsPCASNPLoading(pc, RPgenofile, num.thread=snakemake@threads, verbose=TRUE)
saveRDS(snpL, file = snakemake@output[["Loads"]])

snpgdsClose(RPgenofile)




          
