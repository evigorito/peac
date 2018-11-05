
## prepares a list of files to run fastqc, the first part of this file is a backup for the script with the same name in emedlab, then the script is run in hpc


##########################################################
########################### HPC ##########################
##########################################################

source('/home/ev250/emedlab/snake-pipe/Rfunctions/inputs.fastqc.R')
library(gdata)
library(xtable)
##input:

samp <- fread('/mrc-bsu/scratch/ev250/emedlab/info/PEAC_sampleinfo_path2fq.txt')
samp[!is.na(genotype_dir), GenotypeDirHpc:="/mrc-bsu/scratch/ev250/emedlab/DNA/imputed"]

## some blood samples have bSampleID..QMUL.ID.only but SampleID..QMUL.ID.only is in the name, without "b" and the SampleID..QMUL.or.Genentech. starts with SAM... For snakemake pipeline I need a unique ID in the file name. To disambiguete file names with same SampleID..QMUL.ID.only I replace in the file name the SampleID..QMUL.ID.only.

## Ambry samples have sometimes truncated SampleID..QMUL.or.Genentech. ID in the name. I need to have full name, however for one of the Ambry samples, QMUL2011071, which has 6 files, this will create 2 duplicated file names. I disambiguate by adding "_1" when file name already exists.

## I write the new basenames with path to a location to create sym.links in new col hpc.sym and create the symlinks

samp <- samp.sym(samp, sym="/home/ev250/emedlab/snake-pipe/fastq")


## running all RA samples, use SampleID..QMUL.ID.only. as individual ID after removing starting "b", instead of HospitalNumber as vcf files are named with QMUL identifier, place QMUL id in second col
samp.s <-sel.sam(samp,diagnosis="RA")

## exclude samples SAM9185509, QMUL2010054, SAM9103834 because failed alignment with STAR.

w <- which(samp.s$SampleID..QMUL.or.Genentech. %in% c('SAM9185509', 'QMUL2010054', 'SAM9103834'))

## save as csv to avoid clashes with "9mth (3mths late)"
write.table(samp.s[-w,],"/mrc-bsu/scratch/ev250/emedlab/RNA/objects/RA.csv", row.names=F, quote=F, col.names=F, sep=",")


## ID for GT is QMUL for fastq files for  from synovium or blood sample baseline, probably to record that GT was done from synovium or blood sample at baseline.

###################################################################################
## Adding ethnicity, gender and GT identifier to all samples
###################################################################################

## 141 samples
g.et <- fread('/mrc-bsu/scratch/ev250/emedlab/info/PEAC_genotype_gender_ethnicity.csv')

##g.et2 <- read.xls('/mrc-bsu/scratch/ev250/emedlab/info/SamplePopulations.xlsx')

samp.w <- sel.sam(samp)

samp.et <- merge(samp.w,g.et[,1:3,with=F], by.x="SampleID..QMUL.ID.only.", by.y="qmul_id", all.x=T)

## get the sample names as in vcf to link GT to individual

## bash: bcftools query -l PEAC_1000G_Phase3_Oct14_chr1.vcf.gz > sample_id.txt

vcfID <- fread('/mrc-bsu/scratch/ev250/emedlab/DNA/imputed/sample_id.txt',header=F)

## 138 samples in vcf file, "QMUL2011018" "QMUL2012057" "QMUL2013074" missing from g.et

## only 93 in my dataset

## merge vcf_id with samp.et

## need to create a column that removes starting "b" from SampleID..QMUL.ID.only. to match vcf_id

samp.et[, vcf_id:=gsub("^b", "", SampleID..QMUL.ID.only.)]

samp.et  <-  merge(samp.et, vcfID, by.x="vcf_id", by.y="V1", all.x=T)

## I remove from vcf_id those entries which are not in vcf_id, only keep the ones I have genotype information in the vcf

samp.et[ ,vcf_id:=ifelse(vcf_id %in% vcfID$V1, vcf_id, NA)]

setkey(samp.et,HospitalNumber)

## individual ID is hospital number. Make sure to add  Ethnicity, Gender and vcf_id in all samples from the same individual

## check if missing values for vcf_id are missing for ethnicity, see next step

sum(is.na(samp.et[!is.na(vcf_id), Ethnicity])) == 0
[1] TRUE

## get unique pairs of ID and ethnicity (vcf_id)

pairs <- unique(samp.et[!is.na(HospitalNumber) & !is.na(Ethnicity),.(HospitalNumber,Gender,  Ethnicity, vcf_id)])

setkey(samp.et,HospitalNumber,Gender,  Ethnicity, vcf_id)

for(i in 1:nrow(pairs)){
    samp.et[HospitalNumber==pairs$HospitalNumber[i] & is.na(Ethnicity), c('Gender',  'Ethnicity', 'vcf_id') := lapply(2:4, function(j) unlist(pairs[i,j,with=F]))]
}

## Simplify GT info, add column for yes or no for GT info

samp.et[, QCd.GT:= ifelse(is.na(GenotypeDirHpc), "NO", "YES")]

## save samp.et
write.table(samp.et, "/mrc-bsu/scratch/ev250/emedlab/RNA/objects/PEAC_eth.txt", row.names=F)





################################
## Questions for QMUL
###############################

## Starting point: look at eQTL at Baseline for RA synovium or blood, using Genentech batches.

## HOW TO GROUP SAMPLES?? Ideally when working with noGT for calling variants from RNA we want to pool as much high quality information per individual (tissue, timepoint) as possible.

# How to deal with paired vs single, just focus on paired?

## For calling variants I need for each fastq sample:
<individual_id><sample_name><sequencer_id><flowcell_id><lane_number><adpater_seq>

## excluded samples SAM9185509, QMUL2010054, SAM9103834 because failed alignment with STAR.
    
## Sample names in vcf files is QMUL_id, I am linking to HospitalNumber to determine whether samples come from same individual.

## For running our method with no genotype data we strongly rely on ethnicity. I found some data that seems to be self-reported in file PEAC_genotype_gender_ethnicity.csv'

## missing ethnicity by diagnosis
samp.et[is.na(Ethnicity),.N, Diagnosis]
                      Diagnosis  N
1:                           NA 55
 2:                           UA 37
 3:                           RA 28
 4:                          PsA 28
 5:                         Mono  7
 6:                        UA/RA  1
 7:                       UAPsA?  1
 8:                          TBC  2
 9:       Spondyloarthritis(SPA)  1
10: CTD(connectivetissuedisease)  2
11:                SLE(wthdrawn)  1
12:                       RAPsA?  2


exp.mis <-'is.na(Diagnosis) & Batch %in% grep("Genen",unique(Batch),value=T) '
vars <- c('Tissue', 'Reads', 'Batch', 'QCd.GT')
ord <- c(vars, "N")
mis.samp <- samp.et[eval(parse(text=exp.mis)), .N, vars]

print(xtable(mis.samp), include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")

## The 55 samples with Diagnosis NA from above were run in GenentechBatch2/3 Paired reads from blood or synovium BUT with NA timpoint, HospitalNumber,  Gender. 
samp.et[is.na(HospitalNumber),]

## sample QMUL2011071, was the sample split in 6 or the last 2 are repeating the first 2???
 SampleID..QMUL.ID.only. SampleID..QMUL.or.Genentech. HospitalNumber Batch
1:             QMUL2011071                  QMUL2011071       WH244100 Ambry
    Reads Timepoint   Tissue Diagnosis GenotypeDirHpc
1: Paired  Baseline Synovium        UA             NA
                                                                 hpc.sym.1
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_TAGGC_L003_R1.fastq.gz
                                                                 hpc.sym.2
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_TAGGC_L003_R2.fastq.gz
                                                                 hpc.sym.3
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_TAGGC_L007_R1.fastq.gz
                                                                 hpc.sym.4
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_TAGGC_L007_R2.fastq.gz
                                                                   hpc.sym.5
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_1_TAGGC_L003_R1.fastq.gz
                                                                   hpc.sym.6
1: /home/ev250/emedlab/snake-pipe/fastq/QMUL2011071_1_TAGGC_L003_R2.fastq.gz
   Gender Ethnicity peac_id
1:     NA        NA      NA

################################ RA samples ###################################

## excluded samples SAM9185509, QMUL2010054, SAM9103834 because failed alignment with STAR.

w <- which(samp.et$SampleID..QMUL.or.Genentech. %in% c('SAM9185509', 'QMUL2010054', 'SAM9103834'))
exc <- samp.et[w, c(1:8,17)]

## split DT into various tables when too many columns 
# define a function that takes two parameters:https://stackoverflow.com/questions/31289987/r-xtable-wrap-overflowing-columns-into-subtables
# - your long data.frame 
# - the number of columns you want to print in one table
varTable <- function( agg, cols ) 
{
  tables <- ceiling( length( agg ) / cols )    # number of tables to produce
  # list, first element is number of sub-tables
  p <- list( tables )
  # produce as many table as needed with the full number of columns
  for( i in 0 : ( tables - 2 ) ) p[[ i + 2 ]] <- xtable( agg[, ( cols * i + 1):( cols * i + cols ), with=F ] )
  # last table may have less columns and takes the caption
  p[[ i + 3 ]] <- xtable( agg[ , ( cols * ( i + 1  ) + 1):( length( agg ) ), with=F ]  )
  # return the list with xtable objects that can now be printed one by one
  return( p )
}

var <- varTable(agg=exc,cols= 3)

for( i in 2 : ( var[[ 1 ]] + 1 ) ) print( var[[ i ]], hline.after = 0 ,include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")





##################

samp.et <- samp.et[-w,]

## get number of samples and number of individuals per groups
## are blood paired also single for same ind? NO


exp <- 'Diagnosis=="RA" & Batch %in% grep("Genen",unique(Batch),value=T) & Timepoint == "Baseline" '
vars <- c('Tissue', 'Reads', 'QCd.GT')
ord <- c(vars, "N")
t1 <- samp.et[eval(parse(text=exp)), .N, vars][order(get(ord)),]
print(xtable(t1), include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")

t2 <- samp.et[eval(parse(text=exp)),][Tissue=="Synovium", .N, by=Ethnicity]

t.split <- Reduce(cbind, lapply(seq(1,24,6), function(i) t2[i:(i+5),]))
print(xtable(t.split, NA.string="NA") , include.rownames = FALSE,  booktabs = TRUE, size="scriptsize", NA.string="NA")


exp <- 'Diagnosis=="RA" & Batch %in% grep("Genen",unique(Batch),value=T) '
vars <- c("Tissue","Reads","Timepoint", "Ethnicity", "HospitalNumber")
## samples
samp.et[eval(parse(text=exp)), .N, by=.(Tissue,Reads,Timepoint, is.na(Ethnicity))][order(Tissue,Reads,Timepoint,N),]

      Tissue  Reads         Timepoint is.na  N
 1:    Blood Paired              6mth FALSE  9
 2:    Blood Paired          Baseline  TRUE  6
 3:    Blood Paired          Baseline FALSE 52
 4:    Blood Single              6mth  TRUE  2
 5:    Blood Single              6mth FALSE 42
 6:    Blood Single          Baseline FALSE 10
 7: Synovium Paired              6mth  TRUE  6
 8: Synovium Paired              6mth FALSE 48
 9: Synovium Paired 9mth (3mths late)  TRUE  1
10: Synovium Paired          Baseline  TRUE  8
11: Synovium Paired          Baseline FALSE 86

vars2 <- c("Tissue","Reads","Timepoint", "QCd.GT", "HospitalNumber")
samp.et[eval(parse(text=exp)), .N, by=eval(vars2[1:4])][order(Tissue,Reads,Timepoint,N),]

print(xtable(samp.et[eval(parse(text=exp)), .N, by=eval(vars2[1:4])][order(Tissue,Reads,Timepoint,N),]), include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")


## Duplicated RA Genentech RNA-seq samples
dup <- duplicated(samp.et[eval(parse(text=exp)),vars,with=F])
if(any(dup)) {
    dup <- samp.et[eval(parse(text=exp)),][dup,HospitalNumber]
    samp.et[HospitalNumber %in% dup ,c("SampleID..QMUL.or.Genentech.",  vars), with=F]
}

  SampleID..QMUL.or.Genentech.   Tissue  Reads Timepoint  Ethnicity
1:                  SAM20389187 Synovium Paired      6mth Bangladesh
2:                  SAM24297997 Synovium Paired      6mth Bangladesh   
   HospitalNumber
1:      WH1326872
2:      WH1326872

dup.t <- samp.et[HospitalNumber %in% dup ,c( "SampleID..QMUL.or.Genentech.", "SampleID..QMUL.ID.only.","Batch", vars2), with=F][1:2,]

var.dup <- varTable(agg=dup.t,cols= 3)

for( i in 2 : ( var.dup[[ 1 ]] + 1 ) ) print( var.dup[[ i ]], hline.after = 0 ,include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")

## Ambry samples

amb.q <- samp[SampleID..QMUL.or.Genentech. %in% c('QMUL2010009', 'QMUL2011071'),c(2:5,8:9,10), with=F]
var.amb <- varTable(agg=amb.q,cols= 6)

for( i in 2 : ( var.amb[[ 1 ]] + 1 ) ) print( var.amb[[ i ]], hline.after = 0 ,include.rownames = FALSE,  booktabs = TRUE, size="scriptsize")


exp <- 'Diagnosis=="RA" & !Batch %in% grep("Genen",unique(Batch),value=T)'
samp.et[eval(parse(text=exp)), .N, by=.(Tissue,Reads,Timepoint, is.na(Ethnicity))][order(Tissue,Reads,Timepoint,N),]

     Tissue  Reads Timepoint is.na  N
1: Synovium Paired      6mth FALSE  2
2: Synovium Paired  Baseline  TRUE  5
3: Synovium Paired  Baseline FALSE 18

dup <- duplicated(samp.et[eval(parse(text=exp)),vars,with=F])
if(any(dup)) {
    dup <- samp.et[eval(parse(text=exp)),][dup,HospitalNumber]
    samp.et[HospitalNumber %in% dup ,c("SampleID..QMUL.or.Genentech.", vars), with=F]
}

##################################################### No RA ############################################################
exp <- 'Diagnosis!="RA" & Batch %in% grep("Genen",unique(Batch),value=T)'

## samples Genentech
samp.et[eval(parse(text=exp)), .N, by=.(Tissue,Reads,Timepoint, is.na(Ethnicity))][order(Tissue,Reads,Timepoint,N),]
     Tissue  Reads Timepoint is.na  N
1:    Blood Paired      6mth  TRUE  5
2:    Blood Paired  Baseline  TRUE 16
3: Synovium Paired      6mth  TRUE 17
4: Synovium Paired  Baseline FALSE  1
5: Synovium Paired  Baseline  TRUE 38

## Duplicated samples?
dup <- duplicated(samp.et[eval(parse(text=exp)),c("Diagnosis",vars),with=F])
if(any(dup)) {
    dup <- samp.et[eval(parse(text=exp)),][dup,HospitalNumber]
    samp.et[HospitalNumber %in% dup ,c("SampleID..QMUL.or.Genentech.", "Diagnosis", vars), with=F]
}

   SampleID..QMUL.or.Genentech.   Tissue  Reads Timepoint Ethnicity
1:                  SAM24298044 Synovium Paired  Baseline        NA
2:                  SAM24298084 Synovium Paired  Baseline        NA
   HospitalNumber
1:        6479966
2:        6479966



exp <- 'Diagnosis!="RA" & !Batch %in% grep("Genen",unique(Batch),value=T)'

## samples Ambry
    Tissue  Reads Timepoint is.na N
1: Synovium Paired  Baseline  TRUE 6

## Duplicated samples? NO
dup <- duplicated(samp.et[eval(parse(text=exp)),c("Diagnosis",vars),with=F])
if(any(dup)) {
    dup <- samp.et[eval(parse(text=exp)),][dup,HospitalNumber]
    samp.et[HospitalNumber %in% dup ,c("SampleID..QMUL.or.Genentech.", "Diagnosis", vars), with=F]
}

#############################################################################

## Look at ethnicity diversity, one example for RA
samp.et[Diagnosis=="RA" & Batch %in% grep("Genen",unique(Batch),value=T) & Reads=="Paired" & Timepoint=="Baseline", .N, by=.(Tissue, Ethnicity)][order(Tissue,N),]
      Tissue       Ethnicity  N
 1:    Blood          Indian  1
 2:    Blood      Afro-Carib  1
 3:    Blood      white brit  1
 4:    Blood Asian-Pakistani  1
 5:    Blood       Caribbean  1
 6:    Blood           asian  1
 7:    Blood      Mixed Brit  1
 8:    Blood         Chinese  1
 9:    Blood       Caucasian  2
10:    Blood      Black Brit  2
11:    Blood           White  2
12:    Blood      White Euro  2
13:    Blood           white  2
14:    Blood      Black Afri  3
15:    Blood           Asian  4
16:    Blood           Black  5
17:    Blood              NA  6
18:    Blood      Bangladesh  7
19:    Blood      White Brit 15
20: Synovium        Sudanese  1
21: Synovium      white brit  1
22: Synovium      Black Cari  1
23: Synovium Asian-Pakistani  1
24: Synovium         African  1
25: Synovium     White other  1
26: Synovium       Caribbean  1
27: Synovium        Filipino  1
28: Synovium      Asian-Japa  1
29: Synovium      Mixed Brit  1
30: Synovium           White  2
31: Synovium      White Euro  2
32: Synovium           white  2
33: Synovium       Pakistani  2
34: Synovium          Indian  3
35: Synovium      Afro-Carib  3
36: Synovium      Black Afri  4
37: Synovium      Black Brit  5
38: Synovium       Caucasian  6
39: Synovium           Black  6
40: Synovium      Bangladesh  7
41: Synovium              NA  8
42: Synovium           Asian 10
43: Synovium      White Brit 24







###################################################################################
################################ eMedLab #########################################

## add full path to fasq file to meta.data of PEAC samples

source('/home/elenav/projects/eqtl/snake-pipe/Rfunctions/inputs.fastqc.R')

x <-  '/mnt/volume/PEAC/RNAseq/RawData/'
y <- '/mnt/volume/PEAC/RNAseq/PEAC_sampleinfo.csv'
z <- '/mrc-bsu/scratch/ev250/emedlab/RNA/fastq/'

## some samples from GenentechBatch1 are in /mnt/volumme1

samp <- meta.fq(x,y,z,batch=NULL, w='/mnt/volume1/')

## add column counting number of fq files per sample (entries for same id)
l <- samp[,.N, by=.(SampleID..QMUL.or.Genentech.)]
for(i in unique(samp$SampleID..QMUL.or.Genentech.)){
    samp[SampleID..QMUL.or.Genentech.==i, N.files:=seq(1,l[SampleID..QMUL.or.Genentech.==i,N])]
}


## add column with path to dir with vcf genotype files per chromosome

geno.samples <- as.vector(names(fread('/home/elenav/projects/eqtl/snake-pipe/objects/PEAC.samples.genotype.txt')))

samp[gsub("b", "", SampleID..QMUL.ID.only.) %in% geno.samples, genotype_dir:='/mnt/volume/PEAC/RNAseq/Genotype/imputed_vcf']

## individual ID is hospital number. Make sure to add genotype_id in all samples from the same individual

## get unique pairs of ID and genotype_id

pairs <- unique(samp[!is.na(HospitalNumber) & !is.na(genotype_dir),.(HospitalNumber,genotype_dir)])

setkey(samp,HospitalNumber,genotype_dir)

for(i in 1:nrow(pairs)){
    samp[HospitalNumber==pairs$HospitalNumber[i], genotype_dir:=pairs$genotype_dir[i]]
}


## save
write.table(samp,file='/home/elenav/projects/eqtl/snake-pipe/objects/PEAC_sampleinfo_path2fq.txt', row.names=F)
## save
write.table(samp,file='/home/elenav/projects/eqtl/snake-pipe/objects/PEAC_sampleinfo_path2fq.txt', row.names=F)

## compare files found in /mnt or /data, softlink to home
a <- fread('/home/elenav/projects/eqtl/snake-pipe/objects/PEAC_sampleinfo.csv')
b <- fread('/home/elenav/projects/eqtl/snake-pipe/objects/PEACsampleinfo.csv')
c <- fread('/home/elenav/projects/eqtl/snake-pipe/objects/metadata.csv')




