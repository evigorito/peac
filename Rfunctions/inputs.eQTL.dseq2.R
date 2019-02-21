source('/home/ev250/Bayesian_inf/trecase/Functions/real.data.R')

##library(devtools)

##install_github("chr1swallace/GUESSFM", ref="groups")
library(GUESSFM)


#' Prepare inputs for eQTL with DEseq2
#'
#' This function allows you to run Btrecase for one gene and multiple pre-selected snps. When there is no enough information to ASE counts or rSNP is not in the reference panel, the function will run bayesian negative binomial model only.
#' @param gene gene id for the gene to run
#' @param chr chromosome where the gene is, example chr=22
#' @param snps either cis-window or character vector with pos:ref:alt allele for each snp, defaults to cis-window
#' @param counts.f path to file with filtered counts: rows genes, first col gene_id followed by samples, prepared in inputs.R
#' @param gene.coord path to file listing gene coordinates and exons, prepared in input.R
#' @param vcf path to vcf file with  GT for the chromosome where the gene is
#' @param nhets minimun number of het individuals in order to run the minumn model (NB only), defaults to 5
#' @param tag.threshold numeric with r2 threshold (0-1) for grouping snps to reduce the number of running tests, to disable use "no"
#' @param out path to save outputs, default to current directory
#' @param prefix optional prefix for saving tables, if NULL gene_id.eqtl will be used
#' @param missing numeric whith maximun percentage of missing values allowed, defaults to 5%
#' @keywords inputs negative binomial known genotype regulatory DESEQ2
#' @export
#' @return saves a list with inputs for running eQTL with DESeq2, 2 elements: 1) data table with counts for gene; 2) data table with sample genotypes coded as 0,1,2 and NA plus columns SNP information. Also saves inputs into "out" directory. Saves file to match tags to SNPs, when tag option is chosen.
#' in.deseq2()

in.deseq2 <- function(gene, chr, snps=5*10^5,counts.f,gene.coord,vcf, nhets=5,tag.threshold=.9, out=".", prefix=NULL, missing=5) {
    
    ## check inputs
    files <- c(counts.f,gene.coord,vcf)
    w <- !file.exists(files)
    if(any(w)) stop(paste("invalid file names ", paste(files[w], collapse= ", ")))
    if(!(tag.threshold >=0 & tag.threshold <=1 | tag.threshold =="no")) stop("invalid  'tag.threshold' argument")
    if(!(dir.exists(out))) stop("invalid 'out' argument")
    num <- list(nhets, missing)
    names(num) <-  c("nhets", "missing")
    w <- !sapply(num, is.numeric)
    if(any(w)) stop(paste("invalid arguments:",paste(names(num)[w], collapse=", ")))
    
    ## Extract inputs for gene

    ## get counts 
    counts.g <- fread(cmd=paste("grep -e gene_id -e ",gene,counts.f), header=TRUE)
    if(nrow(counts.g)==0) stop("Gene id is not found in count matrix")
    counts.g <- counts.g[,2:ncol(counts.g),with=F] ## removes gene_id
    
    ## get rsnps and extract GT, remove non-informative snps (missing or hom in all samples)
    ## create dt to collect rsnps excluded from analysis
    rsnps.ex <- data.table(id=character(), reason=character())
    gcoord <- fread(gene.coord)

    
    if(is.numeric(snps)) {
        cis_window <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + c(-snps, snps)},
                               error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        
        gt.as <- vcf_w2(vcf,chr, st=cis_window[["start"]], end=cis_window[["end"]], f.arg='"%CHROM %POS %ID %REF %ALT[ %GT]\\n"', exclude="yes")
        if(is.character(gt.as)) stop(print(gt.as))
        rsnps.ex <- gt.as$excluded
        gt.as <- gt.as$keep
        rs <- copy(gt.as)
        
    } else {
        pos <- as.numeric(sapply(strsplit(snps, split=":"), `[[`,1))
        w <- which(!is.na(pos))
        if(!length(w)) stop(cat("Invalid format for snps ", snps[w]))
        ## get gene start and end, ciswindow=0
        st_end <- tryCatch({gcoord[gene_id==gene & chrom==chr,.(start,end)] + rep(0, 2)},
                               error=function(e) {paste("Gene " ,gene, "and chromosome", chr, "are incompatibles in gene.coord input")})
        if(is.character(st_end)) stop(st_end)
        ## construct cis-window with snps, making sure to include the whole gene
        m <- min(pos) < st_end[[1]]
        M <- max(pos) > st_end[[2]]
        cis_window <- ifelse(m, min(pos), st_end[[1]])
        cis_window <- c(cis_window, ifelse(M, max(pos), st_end[[2]]))
        gt.as <- vcf_w2(vcf,chr, cis_window[1], cis_window[2], f.arg='"%CHROM %POS %ID %REF %ALT[ %GT]\\n"', exclude = "yes")
        if(is.character(gt.as)) stop(print("snps not found in vcf"))
        rsnps.ex <- gt.as$excluded[id %in% snps,]
        gt.as <- gt.as$keep
        rs <- gt.as[id %in% snps,]
        if(!nrow(rs)) stop("Missing GT or homozygous snps in all samples")        
    }    
    ## further process of rsnps   
    ## recode to 0,1,2 scale, missing data NA
    rec.rs <- rec_unphase_rSNPs(y=rs)
    ## remove those with less than nhets
    GT.aux <- rec.rs[,grep("_GT",names(rec.rs)),with=F] ## to make easier calculation of correlations, etc.   
    ## counts number of hets per rsnp
    rec.rs[, nhet:=apply(GT.aux ,1, function(i) sum(i==1, na.rm=T))] 
    ## remove snps with less than min hets
    w <- rec.rs[nhet<nhets, which = TRUE]
    rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[w,id], reason=rep(paste("rsnp with less than", nhets ,"het ind."), length(w))))
    rec.rs <- rec.rs[!w,]
    if(nrow(rec.rs)==0) stop(cat("No rsnp with at least", nhets ,"het ind."))
    ## add number of missing GT per snp
    rec.rs[, mis.GT:= apply(GT.aux[!w,] ,1, function(i) sum(is.na(i)))]
    
    ## remove SNPS with more missing values than allowed
    rem <- rec.rs[mis.GT>=ncol(GT.aux)*missing/100, which=TRUE]
    rsnps.ex <- rbind(rsnps.ex, data.table(id=rec.rs[rem,id], reason=rep(paste("rsnp with more missing genotypes than", missing ,"%."), length(rem))))
    rec.rs <- rec.rs[!rem,]
    
    if(tag.threshold!="no") {
        ## Group rsnps by r2, recode rec.rs for input in tags function from GUESSFM
        if(nrow(rs)==1) stop("Only one regulatory snp to test, please set tag.threshold='no' \n Cannot cluster one snp only")
        re.guess <- rec.guess(DT=rec.rs)
        x <- as(re.guess-1, "SnpMatrix")
        rtag <- GUESSFM::tag(X=x,tag.threshold=tag.threshold)
        ## save rtag as data.table 
        dt <- data.table(Gene_id=gene,tag=tags(rtag), SNP=rtag@.Data)
        if(!is.null(prefix)){
            write.table(dt,paste0(out,"/",prefix,".tags.lookup.txt"), row.names=FALSE)
        } else {
            write.table(dt,paste0(out,"/",gene,".eqtl.tags.lookup.txt"), row.names=FALSE)
        }
        
        ## restrict rsnp to tag snps
        rec.rs <- rec.rs[id %in% unique(tags(rtag)),]
    }
    ## save inputs
    l <- list(counts=counts.g, genotype=rec.rs)
    saveRDS(l, paste0(out, "/",gene,".dseq2.inputs.rds")

}

