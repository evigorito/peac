
## prepares a list of files to run fastqc, the first part of this file is a backup for the script with the same name in emedlab, then the script is run in hpc



########################### HPC ##########################
source('/home/ev250/emedlab/snake-pipe/Rfunctions/inputs.fastqc.R')

##input:

samp <- fread('/mrc-bsu/scratch/ev250/emedlab/info/PEAC_sampleinfo_path2fq.txt')
samp[!is.na(genotype_dir), GenotypeDirHpc:="/mrc-bsu/scratch/ev250/emedlab/DNA/imputed"]

## some blood samples have bSampleID..QMUL.ID.only but SampleID..QMUL.ID.only is in the name, without "b" and the SampleID..QMUL.or.Genentech. starts with SAM... For snakemake pipeline I need a unique ID in the file name. To disambiguete file names with same SampleID..QMUL.ID.only I replace in the file name the SampleID..QMUL.ID.only.

## Ambry samples have sometimes truncated SampleID..QMUL.or.Genentech. ID in the name. I need to have full name, however for one of the Ambry samples, QMUL2011071, which has 6 files, this will create 2 duplicated file names. I disambiguate by adding "_1" when file name already exists.

## I write the new basenames with path to a location to create sym.links in new col hpc.sym and create the symlinks

samp <- samp.sym(samp, sym="/home/ev250/emedlab/snake-pipe/fastq")


## running all RA samples
samp.s <-sel.sam(samp,diagnosis="RA")

## exclude samples SAM9185509, QMUL2010054, SAM9103834 because failed alignment with STAR.

w <- which(samp.s$SampleID..QMUL.or.Genentech. %in% c('SAM9185509', 'QMUL2010054', 'SAM9103834'))

## save as csv to avoid clashes with "9mth (3mths late)"
write.table(samp.s[-w,],"/mrc-bsu/scratch/ev250/emedlab/RNA/objects/RA.csv", row.names=F, quote=F, col.names=F, sep=",")

#for testing pipeline
##write.table(samp.s[1:2,],"/mrc-bsu/scratch/ev250/emedlab/RNA/objects/RA.csv", row.names=F, quote=F, col.names=F, sep=",")



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




