library(data.table)
library(GenomicAlignments)
library(GenomicFeatures)

#' prepares exon by gene input file for calculating counts per gene and prepares data table gene_coordinates
#'
#' @param gtf gtf file with genome reference annotations
#' @param out.ebg full path and name to save "GRangesList" with names genes and ranges exons
#' @param out.gene.coord full name to save text file with gene coordinates
#' @return saves a "GRangesList" with names genes and ranges exons and a text file with gene_id start and end, if in plus strand or otherway round if in minus strand.
#' @export
#'
#' gtf_various()
gtf_various<- function(gtf,out.ebg, out.gene.coord){

    txdb <- makeTxDbFromGFF(gtf, format="gtf", circ_seqs=character())
    ebg <- exonsBy(txdb, by="gene")
    saveRDS(ebg, out.ebg)

    gene_id <- names(ebg)
    st <- min(start(ebg))
    end <- max(end(ebg))
    dt <- data.table(gene_id=gene_id, start=st, end=end)
    gene_id <- names(ebg)
    st <- min(start(ebg))
    end <- max(end(ebg))
    dt <- data.table(gene_id=gene_id, start=st, end=end)

    ## add chrom to dt:
    chrom <- select(txdb, keys= gene_id, columns = "TXCHROM", keytype = "GENEID")

    dt <- merge(dt, chrom, by.x="gene_id", by.y="GENEID", sort=F)
    setnames(dt, "TXCHROM", "chrom")
    write.table(dt, file=out.gene.coord, row.names=FALSE)
}

## running function
gtf_various(snakemake@input[[1]], snakemake@output[[1]], snakemake@output[[2]])

