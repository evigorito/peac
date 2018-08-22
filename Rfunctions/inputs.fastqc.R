### Prepare inputs for Snakefile

library(data.table)

################################ hpc ######################################


#' Selects samples from PEAC metadata to run in python pipeline, converts to wide format
#
#' @param x data table with PEAC metadata
#' @param batch character vector with batch name to select, if null selects all
#' @param reads whether Paired or  Single, if null selects all
#' @param timepoint whether to select samples at Baseline, 6mth or 9mth**, null includes all 
#' @param tissue character vector with tissue, Blood,Synovium,HealthySynovium, null includes all
#' @param diagnosis character vector with diagnosis: options are RA,UA,PsA,CTD(connectivetissuedisease),Mono,TBC,RAPsA?,SLE(wthdrawn),RA(withdrawn),UA/RA,UAPsA?,Spondyloarthritis(SPA)
#' @export
#' @return data table with selected variables in wide format
#' sel.sam()

sel.sam <- function(x,batch=NULL,reads=NULL,timepoint=NULL,tissue=NULL,diagnosis=NULL) {
    var <- c("batch","reads","timepoint","tissue","diagnosis")
    f.names <- c("Batch","Reads","Timepoint","Tissue","Diagnosis")
    rem <- sapply(var,function(i) length(which(is.null(get(i))))==0)
    DT <- copy(x)
    if(any(rem)){
        var <- var[rem]
        f.names2 <- f.names[rem]       
        for(i in seq_along(var)){
            DT <- DT[get(f.names2[i]) %in% get(var[i]),]
        }
    } 
    ## convert to wide format, select columns for Python
    DT <- reshape(DT[,.(SampleID..QMUL.or.Genentech.,N.files,hpc.sym)], idvar = "SampleID..QMUL.or.Genentech.", timevar = "N.files", direction = "wide")

    ## Add individual ID (HospitalNumber) columns in f.names and GenotypeDirHpc that may become handy when processing pipeline
    DT <- merge(DT, unique(x[,c("SampleID..QMUL.or.Genentech.", f.names,"HospitalNumber", "GenotypeDirHpc"), with=F]), by="SampleID..QMUL.or.Genentech.")

    ## Reorder columns

    setcolorder(DT, c("SampleID..QMUL.or.Genentech.", "HospitalNumber", f.names, "GenotypeDirHpc", grep("hpc", names(DT), value=T)))
    
    return(DT)

}


#' Disambiguate Ambry file names as in hpc.path and allows to have a unique ID which is in the file name for some of the blood samples of PEAC study. Create sym link to a designated location with the new basenames and write this location to the orginal data frame.
#
#' @param x data table with PEAC metadata
#' @param sym character vector with dir path to create sym links
#' @export
#' @return data table with hpc.sym column with full path to sym link and creates the sym link
#' samp.sym()

samp.sym <- function(x,sym) {
    ##identify files where SampleID..QMUL.or.Genentech. ID is not in file name
    ## x[,id.f:= mapply(function(i,j) { length(grep(gsub('^[A-Z]', '',i),basename(j))) } , samp$SampleID..QMUL.or.Genentech. , samp$hpc.path)]
    ## ## When id in file, use the same basename for sym link
    ## x[id.f==1, hpc.sym:=paste0(sym, "/",basename(hpc.path))]

    ## when id not in file, replace file name,  but  remove lower case from samp$SampleID..QMUL.ID.only.
    x[,rep.f:=mapply(function(i,j,k) { gsub(gsub('^[a-z]', '',i),j,basename(k)) } , samp$SampleID..QMUL.ID.only., samp$SampleID..QMUL.or.Genentech. , samp$hpc.path)]
    
    ## Ambry has truncated SampleID, I need to gsub everything before [T,G,A,C] with SampleID..QMUL.or.Genentech.

    x[Batch=="Ambry", rep.f:=mapply(function(i,j) {gsub(".*_[T,G,A,C]", paste0("\\1" ,i,"_"), j)} , samp[Batch=="Ambry",SampleID..QMUL.or.Genentech.] , samp[Batch=="Ambry", rep.f])]
    x[Batch!="Ambry", rep.f:=rep.f]

    x[Batch=='Ambry' & duplicated(rep.f), rep.f := mapply(function(i,j) paste0(i, "_1", gsub(i, "", j)) , samp[Batch=='Ambry' & duplicated(rep.f), SampleID..QMUL.or.Genentech.] , samp[Batch=='Ambry' & duplicated(rep.f) , rep.f] ) ]

    x[, hpc.sym:=paste0(sym, "/", rep.f)][,rep.f:=NULL]

    ## create sym links, remove NA
    DT <- x[!is.na(bastion.path),]

    ## check if dir exists, otherwise create it

    if(!dir.exists(sym)) { dir.create(sym) }
    
    f <- mapply(file.symlink , DT$hpc.path, DT$hpc.sym)
    
    return(x)

}









################################# eMedLab back-up ##############################

## Prepare inputs for Snakefile

library(data.table)


#' Creates symbolic links to fastq files to use so all are in one dir to be called by snakefile
#
#' Given the paths to fastq files and sample characteristics, select files that meet criteria and crea
#' te sym links so all are in one place
#' @param x path to dir fastq files, fastq files may be recursively in subfolders depending on batch.
#' @param y path to sample file, when sample characts are used, defaults to NULL
#' @param z path to dir to create symlink
#' @param d optional character vector with diagnosis (disease), defaults to NULL (no selection) 
#' @param batch character vector with batch to look for: Ambry,"GenentechBatch1" "GenentechBatch2" "GenentechBatch3", defualts to NULL which includes all batches
#' @export
#' @return creates symbolic links to specified files in specified directory
#' sel.fq()

sel.fq <- function(x,y,z,d=NULL,batch=NULL){
    ##samples
    samp <- fread(y)
    ## select samples with d
    if(!is.null(d)){
        samp <- samp[Diagnosis==d,]
    }
    
    if(is.null(batch) | "Ambry" %in% batch){

        amb <- ambry.sel.fq(samp)
        tmp <- matrix(unname(c(unlist(amb))), ncol=1)
        ## add column with symlink name 
        tmp <- cbind(tmp, apply(tmp,2,function(i) paste0(z,basename(i))))
        if(length(batch)==1){
            ## create symlik
            f <- lapply(1:nrow(tmp),function(i) file.symlink(tmp[i,1],tmp[i,2]))
            ######################## working here, symlinks need abs path and I am using relative!!!! ###########
            
        }
    }
    if(is.null(batch)){
        rest <- grep("GenentechBatch", unique(samp$Batch), value=T)
    } else {
        rest <- grep("GenentechBatch", batch, value=T)
    }
    fq.rest <- lapply(rest, function(i) genen.sel.fq(samp,batch=i))
    tmp2 <- matrix(unname(c(unlist(fq.rest))), ncol=1)
    ## add column with symlink name 
    tmp2 <- cbind(tmp2, apply(tmp2,2,function(i) paste0(z,basename(i))))
    if(is.null(batch) | "Ambry" %in% batch) {
        tmp3 <- rbind(tmp,tmp2)
        f <- lapply(1:nrow(tmp3),function(i) file.symlink(tmp3[i,1],tmp3[i,2]))
    }
}
    
    


#' Prepare inputs for fastqc by snakefile, sub function for genentech batches
#
#' Deal with Genentech batches in sel.qc
#' @param samp 
#' @param batch character vector with 1 batch to look for: either "GenentechBatch1" "GenentechBatch2" "GenentechBatch3"
#' @export
#' @return list with names of fasq files
#' genen.sel.fq()

genen.sel.fq <- function(samp,batch){
    samp.gen <- samp[Batch==batch,]

    ## match id with fasq.files
    if(batch=="GenentechBatch1" | batch=="GenentechBatch3") {
        fq <- sapply(1:nrow(samp.gen), function(i) list.files(path=paste0(x,samp.gen[i,Batch]), pattern=unique(samp.gen[i,SampleID..QMUL.or.Genentech.]), full.names=TRUE,recursive=TRUE,include.dirs=TRUE))
    
    names(fq) <- samp.gen$SampleID..QMUL.or.Genentech.

    ## get mising samples

    mis.gen <- names(fq[sapply(fq,function(i) !length(i))])

    samp.mis <- samp[SampleID..QMUL.or.Genentech. %in% mis.gen,]

    ## get hits:
    fq <- fq[sapply(fq,function(i) as.logical(length(i)))]
    
    ## All blood samples for GenentechBatch1 are missing
    }

    ## For Genen'Batch2 id is sampleID..QMUL.ID.only. removing the initial b when present (I guess indicates blood samples).
    if(batch=="GenentechBatch2"){

        fq.rescue <- sub("b","",samp.gen$SampleID..QMUL.ID.only.)
        fq.rs <- sapply(fq.rescue, function(i) list.files(path=paste0(x,samp.gen[SampleID..QMUL.ID.only.== i,Batch]), pattern=i, full.names=TRUE,recursive=TRUE, include.dirs=TRUE))
        names(fq.rs) <- samp.gen$SampleID..QMUL.ID.only.

        ## missing samples
        res.mis <- names(fq.rs[sapply(fq.rs,function(i) !length(i))])

        ## get no missing
        fq <- fq.rs[sapply(fq.rs,function(i) length(i)!=0)]
        
    }
    
    return(fq)
}



#' Prepare inputs for fastqc by snakefile, sub function for genentech batches
#
#' Deal with Ambry batch in sel.qc
#' @param samp 
#' @export
#' @return list with names of fasq files
#' ambry.sel.fq()

ambry.sel.fq <- function(samp){

    ## Dir with fasq files is organised differently for Ambry vs rest of samples

    ## Ambry: all samples are in dirs named with sample$id but in some dirs the fasq.gz files do not contain the id name. I search for the dir and within the dir any fasq.gz sample

    samp.a <- samp[Batch=="Ambry",]

    ## dirs with sample name for Ambry are not exactly the same as in id.a, they miss the 3rd character from end to start which is a 0. I need to delete this to get the correct dir names

    dir.pat.a <- sapply(unique(samp.a$SampleID..QMUL.or.Genentech.), function(i) {
        tmp <- unlist(strsplit(i, ""))
        tmp <- paste(tmp[-(nchar(i)-2)], collapse="")
        return(tmp)
    })
    
    dirs.a <- sapply(dir.pat.a, function(i) grep(i, list.dirs(x),value=T))

    ## missing samples
    mis.a <- names(dirs.a[sapply(dirs.a,function(i) !length(i))])

    ## remove missing samples
    dirs.a <- dirs.a[sapply(dirs.a,function(i) length(i)!=0)]

    ## get complete file names

    am.f <- lapply(dirs.a, function(i) list.files(path=i, pattern="fastq.gz",full.names=TRUE))

    return(am.f)
    
} 

#' Add file path to fastq files to PEAC meta data 
#
#' Given the paths to fastq files add them as a column into PEAC meta.data
#' @param x path to dir fastq files in bastion, fastq files may be recursively in subfolders depending on batch.
#' @param y path to sample file, when sample characts are used, defaults to NULL
#' @param z path to dir with fastq files in hpc,basename of files is the same as bastion
#' @param batch character vector with batch to look for: Ambry,"GenentechBatch1" "GenentechBatch2" "GenentechBatch3", defualts to NULL which includes all batches
#' @param w second path to dir fastq files in bastion, /mnt/volume1/ includes fastq files in 'Single' and 'Paired' dirs
#' @export
#' @return data table with extra columns with full path to fasq files to bastion and hpc
#' meta.fq()

meta.fq <- function(x,y,z, batch=NULL,w=NULL){
    ##samples
    samp <- fread(y)
    samp[,V1:=NULL]
    samp <- unique(samp)
    
    if(substr(z, nchar(z), nchar(z)) != "/"){
                z <- paste0(z,"/")
    }
    if(substr(x, nchar(x), nchar(x)) == "/"){ ## to avoid "//" when calling list.files
                x <- substr(x,1,nchar(x)-1)
    }
    
    ## Go by batch, start with Ambry
    ## dirs with sample name for Ambry are not exactly the same as in id.a, they miss the 3rd character from end to start which is a 0. I need to delete this to get the correct dir names
    
    if(is.null(batch) | "Ambry" %in% batch){
        samp.a <- samp[Batch=="Ambry",]
        dir.pat.a <- sapply(unique(samp.a$SampleID..QMUL.or.Genentech.), function(i) {
            tmp <- unlist(strsplit(i, ""))
            tmp <- paste(tmp[-(nchar(i)-2)], collapse="")
            return(tmp)
        })
        ## sample QMUL2010009 is labelled as QMUL20109, correct
        dir.pat.a[which(dir.pat.a == "QMUL201009")] <- 'QMUL20109'

        ## get dirs
        dirs.a <- sapply(dir.pat.a, function(i) grep(i, list.dirs(x),value=T))

        ## get complete file names
        am.f <- rbindlist(lapply(seq_along(dirs.a), function(i) data.table(bastion.path=list.files(path=dirs.a[[i]], pattern="fastq.gz",full.names=TRUE), id=names(dirs.a)[[i]])))
        
        ## Add path to samp
        samp <- merge(samp,am.f, by.x="SampleID..QMUL.or.Genentech." , by.y="id", all.x=T)
        if(!is.null(batch)){
            
            samp[,hpc.path:=paste0(z,basename(bastion.path))]
            return(samp)
        }
    }
    if(is.null(batch)){
        rest <- grep("GenentechBatch", unique(samp$Batch), value=T)
    } else {
        rest <- grep("GenentechBatch", batch, value=T)
    }
    fq.rest <- lapply(rest, function(i) genen.meta.fq(x=x,y=samp,batch=i))
    names(fq.rest) <- rest
    if(is.null(batch) | "GenentechBatch2" %in%  batch){
        ## merge GenentechBatch2, uses "SampleID..QMUL.ID.only."
        samp <- merge(samp,fq.rest$GenentechBatch2,by.x="SampleID..QMUL.ID.only." , by.y="id", all.x=T)
        samp[!is.na(bastion.path.y), bastion.path.x:=bastion.path.y]
        setnames(samp,"bastion.path.x","bastion.path")
        samp[,bastion.path.y:=NULL]
        if( "GenentechBatch2" %in%  batch ){
            samp[,hpc.path:=paste0(z, basename(bastion.path))]
            return(samp)
        }
        ## remove Gene'B2 from fq.rest
        fq.rest$GenentechBatch2<-NULL
    }
    fq.d <- rbindlist(fq.rest)
    if(!is.null(w)){
        if(substr(w, nchar(w), nchar(w)) == "/"){ ## to avoid "//" when calling list.files
            w <- substr(w,1,nchar(w)-1)
        }
        fq.w <- genen.v1(x=samp,y=w)
        fq.d <- rbind(fq.d,fq.w)
    }
        
    ## merge Gene'B1 and/or B3 by  SampleID..QMUL.or.Genentech.
    samp <- merge(samp,fq.d, by.x="SampleID..QMUL.or.Genentech." , by.y="id", all.x=T)
    samp[!is.na(bastion.path.y), bastion.path.x:=bastion.path.y]
    setnames(samp,"bastion.path.x","bastion.path")
    samp[,bastion.path.y:=NULL]
    samp[,hpc.path:=paste0(z, basename(bastion.path))]

    ## remove duplicate entries due to SampleID..QMUL.ID.only having 2 entries in SampleID..QMUL.or.Genentech.
    return(unique(samp))
    
}


#' Add file path to fastq files for PEAC data, sub function for genentech batches
#
#' Deal with Genentech batches in meta.fq
#' @param y data table with metadata from PEAC
#' @param x path to PEAC dir with fastq files
#' @param batch character vector with 1 batch to look for: either "GenentechBatch1" "GenentechBatch2" "GenentechBatch3"
#' @export
#' @return data table with path to bastion and id
#' genen.meta.fq()

genen.meta.fq <- function(x,y,batch){
    samp.gen <- y[Batch==batch,]

    ## match id with fasq.files
    if(batch=="GenentechBatch1" | batch=="GenentechBatch3") {
        fq <- lapply(1:nrow(samp.gen), function(i) list.files(path=paste0(x,"/",samp.gen[i,Batch]), pattern=paste0(unique(samp.gen[i,SampleID..QMUL.or.Genentech.]), ".*fastq.gz"), full.names=TRUE,recursive=TRUE,include.dirs=TRUE))

        names(fq) <- samp.gen$SampleID..QMUL.or.Genentech.

        ## remove missing values
        fq <- fq[sapply(fq,function(i) as.logical(length(i)))]

        g.f <- rbindlist(lapply(seq_along(fq), function(i) data.table(bastion.path=fq[[i]], id=names(fq)[[i]])))

        ## ## Add bastion.path to y
        ## y <- merge(y,g.f, by.x="SampleID..QMUL.or.Genentech." , by.y="id", all.x=T)
        ## y[!is.na(bastion.path.y), bastion.path.x:=bastion.path.y]
        ## setnames(y,"bastion.path.x","bastion.path")
        ## y[,bastion.path.y:=NULL]
        ## ##y[,hpc.path:=paste0(z,"/", basename(path.bastion))]
    }

    ## For Genen'Batch2 id is sampleID..QMUL.ID.only. removing the initial b when present (I guess indicates blood samples).
    if(batch=="GenentechBatch2"){

        fq.rescue <- sub("b","",samp.gen$SampleID..QMUL.ID.only.)
        fq.rs <- lapply(fq.rescue, function(i) list.files(path=paste0(x,"/",samp.gen[SampleID..QMUL.ID.only.== i,Batch]), pattern=i, full.names=TRUE,recursive=TRUE, include.dirs=TRUE))
        names(fq.rs) <- samp.gen$SampleID..QMUL.ID.only.

        ## missing samples
        res.mis <- names(fq.rs[sapply(fq.rs,function(i) !length(i))])

        ## get no missing
        fq <- fq.rs[sapply(fq.rs,function(i) length(i)!=0)]
        
        g.f <- rbindlist(lapply(seq_along(fq), function(i) data.table(bastion.path=fq[[i]], id=names(fq)[[i]])))

        ## ## Add bastion.path to y
        ## y <- merge(y,g.f, by.x="SampleID..QMUL.ID.only." , by.y="id", all.x=T)
        ## y[!is.na(bastion.path.y), bastion.path.x:=bastion.path.y]
        ## setnames(y,"bastion.path.x","bastion.path")
        ## y[,bastion.path.y:=NULL]
        
    }
    
    return(g.f)
}

#' Add file path to fastq files for PEAC data, dir /mnt/volume1
#
#' Deal with Genentech batches in meta.fq, this time with blood samples in /mnt/volume1
#' @param x path to PEAC dir with fastq files
#' @param y data table with metadata from PEAC
#' @export
#' @return data table with path to bastion and id
#' genen.v1()

genen.v1 <- function(x,y){

    fq.p <- lapply(1:nrow(x), function(i) list.files(path=paste0(y,'/Paired'), pattern=paste0(unique(x[i,SampleID..QMUL.or.Genentech.]), ".*fastq.gz"), full.names=TRUE,recursive=TRUE,include.dirs=TRUE))

    names(fq.p) <- x$SampleID..QMUL.or.Genentech.

    ## get no missing
    fq.p <- fq.p[sapply(fq.p,function(i) length(i)!=0)]

    g.f <- rbindlist(lapply(seq_along(fq.p), function(i) data.table(bastion.path=fq.p[[i]], id=names(fq.p)[[i]])))

    fq.s <-  lapply(1:nrow(x), function(i) list.files(path=paste0(y,'/Single'), pattern=paste0(unique(x[i,SampleID..QMUL.or.Genentech.]), ".*fastq.gz"), full.names=TRUE,recursive=TRUE,include.dirs=TRUE))
     names(fq.s) <- x$SampleID..QMUL.or.Genentech.

    ## get no missing
    fq.s <- fq.s[sapply(fq.s,function(i) length(i)!=0)]

    g.fs <- rbindlist(lapply(seq_along(fq.s), function(i) data.table(bastion.path=fq.s[[i]], id=names(fq.s)[[i]])))

    return(rbind(g.f,g.fs))
}

