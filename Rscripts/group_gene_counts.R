library(data.table)

## get args from snakemake

in.files <- unlist(snakemake@input)

out.file <- snakemake@output[[1]]

lib.size.file <-  snakemake@output[[2]]

##filter <- as.numeric(snakemake@params[['filter']])

filter =0

## process inputs
lcounts <- lapply(in.files, fread)

counts <- Reduce(function(...) merge(...,by="gene_id"), lcounts)

lib_size= log(colSums(as.matrix(counts[,2:ncol(counts),with=F])))

lib_size <- matrix(lib_size, ncol=1, dimnames=list(names(counts)[2:ncol(counts)], "lib.s"))

## filter counts, remove genes with 0 counts in all samples
counts <- counts[which(rowSums(counts[, 2:ncol(counts), with=F])>filter),]

## save files

write.table(counts, out.file, row.names=F)

saveRDS(lib_size, lib.size.file)

