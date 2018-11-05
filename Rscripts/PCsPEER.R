library(gdsfmt)
library(SNPRelate)
library(data.table)
library(peer)
library(DESeq2)
library(biomaRt)
library(cqn)

#' Compute PCs from peac samples by projecting into the reference panel and computes PEER factors based on DESeq2 variance stabilization or cqn with GC correction. Prepare inputs for MatQTL
#
#' @param gds gds file name with genotype for samples
#' @param snpL rds object with matrix with SNP loadings to project samples
#' @param n number of eigenvectors and PEER to extract, defaults to 10
#' @param ld ld cut-off to thin genotype matrix
#' @param eaf cut-off to thin genotype matrix
#' @param counts file with count data per gene, output from rule group_gene_counts, first column gene id, then samples
#' @param metadata path to file to peac metadata to relate vcf_id with fasq files id (sample names in expression data)
#' @param prefix, character vector with wildcards for covs from snakemake rule PCs_PEER
#' @param gene.coord file with gene name, chrom, start and end
#' @param out named list with output files in this order; file with all required pcs ; file with all required PEER factors by variance stailization; file with all required PEER factors by cqn; covs N files with increasing number of PCs (1 to n), M files with increasing number of PEER factors for cqn; file with genotypes as required by matrixqtl. Snps are filtered by ld and maf to reduce number of tests; expression matrix for cqn; file with gene id, start and end position as required by matrixqtl; file with snp id, chromosome and position. Names as in snakefile rule PCs_PEER
#' 
#' @keywords extract snps gene
#' @export
#' @return saves multiple files: 1) eigen vectors for samples, 
#' pcs.peer()

pcs.peer <- function(gds,snpL, n=10,  ld, eaf, counts, meta, prefix, gene.coord, out) {

    ## get meta data file to relate ids and covariates:
    meta <- fread(meta)
    
    ## get genotypes
    sampgds <- snpgdsOpen(gds)
                      
    ## get eigenvectors
    snpL <- readRDS(snpL)  
    SL <- snpgdsPCASampLoading(snpL, sampgds)  
    ev <- data.table(sample.id = SL$sample.id)
    for(i in 1:n){
        ev[ , paste0("EV",i) := SL$eigenvect[,i] ]
    }

    ## select samples with expression data and make sure both datasets are in the same order
    expr <- fread(counts)
    matExp <- as.matrix(expr[, 2:ncol(expr)])
    rownames(matExp) <- expr$gene_id

    vcf_in <- meta[SampleID..QMUL.or.Genentech. %in% colnames(matExp), .(SampleID..QMUL.or.Genentech., vcf_id, Batch, Gender)]

    ## select and order ev 
    ev  <- ev[sample.id %in% vcf_in$vcf_id, ]
    ev <- ev[order(match(sample.id, vcf_in$vcf_id))]

    
    ## order matExp as in vcf_in
    matExp <- matExp[, vcf_in$SampleID..QMUL.or.Genentech.]

    ## for eqtl only use genes in chrom 1-22
    ## get gene coordinates and subset
    gene.c <- fread(gene.coord)
    gene.c <- gene.c[chrom %in% 1:22,]

    matExp <- matExp[rownames(matExp) %in% gene.c$gene_id, ]

    ## save gene.c as for matrixqtl gene_location file
    setcolorder(gene.c, c("gene_id", "chrom", "start", "end"))
    write.table(gene.c, row.names=F, file=out[['geneLoc']])
                
    ## get cqn normalised expression and  peer factors    
    ## need gene length and GC content from biomart library, use longest transcript

    ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")

    gc <- data.table(getBM(attributes = c("ensembl_gene_id", "transcript_length", "percentage_gene_gc_content"),
                           filters = "ensembl_gene_id",
                           values = rownames(matExp),
                           mart = ensembl))

    setkey(gc, ensembl_gene_id, transcript_length)

    ## get gc and length for the longest transcript
    gc <- gc[, .SD[.N], ensembl_gene_id]
    setnames(gc, names(gc), c("gene_id", "length", "gc"))

    ## cqn normalise
    cqn.peac <- cqn(counts = matExp, x = gc$gc/100, lengths = gc$length,   verbose=T)

    ## get normalised values
    rpkm.peac <- cqn.peac$y + cqn.peac$offset

    ## save normalised values
    write.table(rpkm.peac, file=out[['expressionCqn']])

    ## run model
    model = PEER()
    PEER_setPhenoMean(model,t(rpkm.peac))
    PEER_setNk(model,n)
    PEER_update(model)
    factors = PEER_getX(model)
    colnames(factors) <- paste0("PEER",1:10)

    ## need to transpose ev and factors and then cbind and save
    ev <- t(ev[,2:ncol(ev)])
    factors <- t(factors)
    
    ## I need to save output covars as in rule: 1 pcs with 1:10 peer factors, then 2 pcs with peer 1:10, etc, rows covariate name, cols samples
    indx=0
    for(i in 1:nrow(ev)){
        for(j in 1:nrow(factors)){
            indx= indx + 1 ## to count files 1:100
            tmp <- data.table(rbind(ev[1:i, ,drop=F], factors[1:j, , drop=F]), keep.rownames=T)
            write.table(tmp, row.names=F, col.names=F, file=out[['covs']][indx])
        }
    }

    ## Same with pcs and batch and sex, add pcs 1:10
    Bsex <- vcf_in[,  .(Batch, Gender)]
    ## recode as numeric
    Bsex[Gender == "M", Gender:="0"][Gender == "F", Gender:="1"]
    b <- unique(Bsex$Batch)
    names(b) <- 0:(length(b)-1)
    for (i in seq_along(b)){
        Bsex[Batch == b[i], Batch:= names(b)[i]]
    }
    ## transpose and make numeric
    mBsex <- t(Bsex)
    mBsex <- apply(mBsex, 2, as.numeric)
    rownames(mBsex) <- names(Bsex)
    
    for( i in 1:nrow(ev)){
        tmp <- data.table(rbind(ev[1:i,, drop=F], mBsex), keep.rownames=T)
        write.table(tmp, row.names=F, col.names=F, file=out[['covfix']][i])

    }
          
    ## Genotypes: LD prune
    gprune <- snpgdsLDpruning(sampgds, maf=eaf, ld.threshold=ld)
    SNPprune <- unlist(gprune)

    ## get genotype matrix for selected snps
    gmat <- snpgdsGetGeno(sampgds, snp.id=SNPprune)

    ## add rownames (sample id) and colnames (snps) to gmat for snp.id
    ## get chrom:pos:ref:alt. snp.id in gds is a an inernal id integer
    ## from 1:N
    
    sampId <- read.gdsn(index.gdsn(sampgds, "sample.id"))
    chrom <- readex.gdsn(index.gdsn(sampgds, "snp.chromosome"),sel=SNPprune)
    pos <- readex.gdsn(index.gdsn(sampgds, "snp.position"), sel=SNPprune)
    allele <- readex.gdsn(index.gdsn(sampgds, "snp.allele"), sel=SNPprune)
    ref <- gsub("/.*", "", allele)
    alt <- gsub(".*/", "", allele)

    rownames(gmat) <- sampId
    colnames(gmat) <- paste(chrom, pos,ref,alt,sep=":")

    ## SNPrelate/gds counts the number of "A" alleles, need to change 2 to 0 and 0 to 2
    g <- copy(gmat)
    g[gmat==2] <- 0
    g[gmat==0] <- 2

    ## transpose g, select samples and order as vcf_in and save
    g <- t(g)
    g <- g[, vcf_in$vcf_id]
    write.table(g, file=out[['geno']], quote=F)

    ## make snp location file as per matrixqtl
    snpDT <- data.table(snp=rownames(g), chr=chrom, pos=pos)
    write.table(snpDT, row.names=F, file=out[['snpLoc']])
       
    snpgdsClose(sampgds)
    
}


out=list(covs=snakemake@output[['covars']],
         geneLoc=snakemake@output[['geneLoc']],
         expressionCqn=snakemake@output[['expressionCqn']],
         covfix=snakemake@output[['covfix']],
         geno=snakemake@output[['geno']],
         snpLoc=snakemake@output[['snpLoc']]
         )

pcs.peer(gds=snakemake@input[['peacGds']], snpL=snakemake@input[['Loads']],
         n=snakemake@params[['Nfactors']],
         ld= snakemake@params[['ld']],
         eaf=snakemake@params[['maf']],
         counts=snakemake@input[['expr']],
         meta=snakemake@input[['peacdata']],
         prefix=snakemake@params[['prefix']],
         gene.coord=snakemake@input[['gene_coord']],
         out=out
         )

