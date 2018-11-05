## Snakefile for emedlab


## https://hpc-carpentry.github.io/hpc-python/17-cluster/
## shell.executable("/bin/bash")
## shell.prefix("module load samtools-1.4-gcc-5.4.0-derfxbk; ")
shell.prefix("source ~/.bashrc; ")

configfile: "config.yaml"

def read_samples():
    """
    Function to get names and fastq paths from a sample file specified
    in the configuration. 
    It works for single or paired sequencing with any number of files per sample. Input file is expected to have ID column followed by as many columns as needed to accommodate all fastq files required, use NA for samples with less files.
    For emedlab I have the following order:
     <SampleID..QMUL.or.Genentech.> <unique_individual_id=HospitalNumber><ID in vcf=SampleID..QMUL.ID.only.><Batch><Reads><Timepoint><Tissue><Diagnosis><GenotypeDir><fastq1.R1_path> <fastq1.R2_path> <fastq2.R1_path> <fastq2.R2_path>. This function produces a dictionary of sample_id keys and values(individualID, Batch, Reads, Timepoint, Tissue, Diagnosis, Dir with genotype, and a sub-list of (fq1, fq2, ...,) 
     Input file prepared in /home/ev250/emedlab/snake-pipe/Rscripts/files4fastqc.R
    """
    with open(config['sample_file'], "r") as f:
        samp_dict = {}
        for line in f:
            words = line.strip().split(",")
            words2 = words[9:]
            while 'NA'in words2:
                words2.remove('NA')
            words = words[:9]
            words.append(words2)
            samp_dict[words[0]] = words[1:]
            # while 'NA'in words:
            #     words.remove('NA')
            # ## make dir, first value tells reads = "paired" or "single"
            # if len(words) == 2:
            #     reads="single"
            # else:
            #     reads="paired"
            # samp_dict[words[0]] = [reads, words[8:(len(words))]]
            ## filter dic if paired/single reads
            # if paired in ("yes" ,"y"):
            #     samp_dict = {key:value for key, value in samp_dict.items() if len(value) > 1 }
            # if paired in ("no","n"):
            #     samp_dict = {key:value for key, value in samp_dict.items() if len(value) == 1 }

    return samp_dict

def group_samples():
    """ Function to group samples by batch (Ambry or any of the Genentech), reads, timepoint, diagnosis and genotype info, choosing one file per individual randomly.
     """
    samp_dict = read_samples()
    ## make list with <Batch><Reads><Timepoint><Tissue><Diagnosis><GenotypeDir> from samp_dict
    sub_values = [v[2:8] for k, v in samp_dict.items()]

    ## get unique combinations of variables
    u_val = [list(x) for x in set(tuple(x) for x in sub_values)]
    ## Prepare dic with key: Batch_reads_tiempoint_tissue_diagnosis_geno, batch Ambry or genentech. Only one sample used if more of the same type available for the same individual
    genen_u_ind = group_helper(samp_dict,u_val,"Genentech")
    ambry_u_ind = group_helper(samp_dict,u_val,"Ambry")
    genen_u_ind.update(ambry_u_ind)
    return(genen_u_ind)


def group_helper(samp_dict, u_val, batch):
    """ Helper function to iterate over Genentech or Ambry. Creates a dictionary with keys Batch_reads_tiempoint_tissue_diagnosis_geno, batch Ambry or genentech, and values sample_name as in fasq
 """
    u_val = [list(x) for x in set(tuple(x[1:]) for x in u_val if x[0].startswith(batch))]
    u_ind = {}
    for x in range(0,len(u_val)):
    ## prepare values
        val = [list(v) for v in set(tuple(v[:8]) for (k,v) in samp_dict.items() if v[3:8]==u_val[x] if v[2].startswith(batch))]
        val = [k for (k,v) in samp_dict.items() for x in val if v[:8] == x ]
        ## prepare key
        key= u_val[x].copy()
        if key[4] != "NA":
            key[4] = "Geno"
        else:
            key.remove("NA")
            
        key.insert(0,batch)
        key = "_".join(key)
        ## remove blank space and conflicting characters
        key=key.replace(" ","").replace("(","").replace(")","")
        u_ind[key] = val
    return(u_ind)

def vcf(path):
    """ Creates dictionary with keys chromosome name and values vcf full name.
    argument is the path a file containing full names to vcf/bcf files
    """
    VcfFiles = open(path).read().splitlines()
    chr = [x.split("chr")[1].split(".")[0] for x in VcfFiles]
    VcfDic = dict(zip(chr, VcfFiles))
    return(VcfDic)      

rule all:
    input:
        ##expand(config['output_dir'] + "/STAR/{read}/{sample}/Aligned.sortedByCoord.out.bam" , zip, sample=read_samples().keys(), read=[item[2] for item in read_samples().values()])
        #expand(config['output_dir'] + "/RNA_counts/groups/{group}.txt" , group=group_samples().keys() ) ,
        #expand(config['output_dir'] + "/RNA_counts/groups/{group}_lib_size.rds" , group=group_samples().keys() )
        #expand(config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz", chrom=vcf(config["ref_bcf"]).keys() ),
        #expand(config['output_dir'] + "/DNA/PEAC_chr{chrom}_sub.vcf.gz", chrom=vcf(config["ref_bcf"]).keys() )
        #expand(config['output_dir'] + "/DNA/RP_chr{chrom}_alt_added.bcf", chrom=merge_by_wildcard().keys() )
        #expand(config['output_dir'] + "/eqtl/neg_binom_gt/{group}_chr{chrom}.txt", group=group_samples().keys() , chrom=vcf(config["geno_vcf"]).keys() )
        # expand(config['output_dir'] + "/matqtl/output/{pcs}.{peerBSex}.txt",
        #        pcs=["pcs" + str(x) for x in range(1,int(config['N factors'])+1)],
        #        peerBSex=["peerCqn" + str(x) for x in range(1,int(config['N factors'])+1)] + ["covSexBatch"])
        expand(config['output_dir'] + "/matqtl/output/pcs0.{peerBSex}.txt",
               peerBSex=["peerCqn" + str(x) for x in range(1,int(config['N factors'])+1)] + ["covSexBatch"])

rule star_index:
    """Create index for alignment using STAR"""
    input:
        config['ref_fasta'],
        config['ref_gtf']
    output:
        config['indices']
    shell:
         "{config[STAR]} "
         " --runThreadN {threads} "
         " --runMode genomeGenerate "
         " --genomeDir {output[0]} "
         " --genomeFastaFiles {input[0]} "
         " --sjdbGTFfile {input[1]} "
         " --sjdbOverhang 100 "  

rule star:
    """ Map paired or sigle end reads using STAR, stores single reads in dir 'single' and paired reads in 'paired' """
    input:
        lambda wildcards: read_samples()[wildcards.sample][7]
    output:
        config['output_dir'] + "/STAR/{read}/{sample}/Aligned.sortedByCoord.out.bam" 
    log:
        "logs/{read}/{sample}.log"
    params:
        index=config['indices'],
        read="zcat"
    threads: 16
    run:
        fq=[input] if isinstance(input, str) else input
        fq1 = ",".join(fq[0:len(fq):2])
        fq2 = ",".join(fq[1:len(fq):2])
        if len(fq2)>0:
            assert len(fq1) == len(fq2), "input-> equal number of files required for paired fasq files"
        input_str =  " ".join([fq1, fq2])
        print(input_str)
        out_dir = [item.replace("Aligned.sortedByCoord.out.bam", "") for item in output]
        shell(
            "{config[STAR]} "
            " --runThreadN {threads} "
            " --genomeDir {params.index} "
            " --readFilesIn {input_str} "
            " --readFilesCommand {params.read} "
            " --outSAMtype BAM SortedByCoordinate "
            " --outFileNamePrefix {out_dir} "
            " --outStd Log "
            " {log}")

rule exon_by_gene:
    """ Get exons per gene as a GRangesList from the gft annotation file used for the alignment, if not already done (to be uses for calculating total raw gene counts). Make a file with gene coordinates to define SNPS within cis-window. """
    input:
        config['ref_gtf']
    output:
        config['ebg'],
        config['output_dir'] + "/gene_inputs/gene_coord.txt"
    script:
         "Rscripts/exon_by_gene.R"


rule total_gene_counts:
    """ Calculate total gene counts from RNA-seq BAM files"""
    input:
        config['ebg'] ,
        lambda wildcards: config['output_dir'] + "/STAR/" + read_samples()[wildcards.sample][2] + "/" + wildcards.sample+ "/Aligned.sortedByCoord.out.bam" 
    params:
       mode="Union",
       ignore_strand="TRUE"
    output:       
        config['output_dir'] + "/RNA_counts/{sample}.txt"
    script:
         "Rscripts/total_gene_counts.R"

rule group_gene_counts:
    """ Group total gene counts by tissue, timepoint, diagnosis and batch and make matrix with log(library size)"""
    input:
        lambda wildcards: expand(config['output_dir'] + "/RNA_counts/{sample}.txt",  sample=group_samples()[wildcards.group])
    params:
        filter=config['filter']
    output:      
        config['output_dir'] + "/RNA_counts/groups/{group}.txt",
        config['output_dir'] + "/RNA_counts/groups/{group}_lib_size.rds"
    script:
         "Rscripts/group_gene_counts.R"

rule vcf_pca:
    """ From vcf input per chromosome removes alleles A/T, T/A, C/G or G/C to avoid switching issues when performing PCA. Because PEAC files for chr1-9 are labelled 01 t0 09"""
    input:
        lambda wildcards: vcf(config["geno_vcf"])[wildcards.chrom]
    output:
        config['output_dir'] + "/DNA/PEAC_chr{chrom}_4PCA.vcf.gz" ,
        config['output_dir'] + "/DNA/PEAC_chr{chrom}_4PCA.vcf.gz.tbi"
    run:
        #for x in wildcards.chrom:
        #    if int(x) <=9 :
        #        inp = vcf(config["geno_vcf"])[x]
        #        tmp = config['output_dir'] + "/DNA/PEAC_chr" + x + "_4PCA.vcf.gz"
         #       shell("bcftools view {inp} | awk '{{gsub(/^0/,\"\"); print}}' | sed 's/##contig=<ID=0/##contig=<ID=/' | bgzip -c {tmp}  " 
          #            "mv {tmp} {inp} ; "
           #           "tabix -p vcf {inp} "
           #     )
        shell(
            "bcftools view {input} -e '"' REF = "A" & ALT = "T" '"'  -Ou | "
            "bcftools view -e  '"'REF = "T" & ALT = "A" '"' -Ou | "
            "bcftools view -e  '"'REF = "C" & ALT = "G" '"' -Ou | "
            "bcftools view  -e '"'REF = "G" & ALT = "C" '"' | "
            "bgzip -c >  {output[0]}; " 
            "tabix -p vcf {output[0]}  "
        )    
    
rule ref_panel_alt:
    """ Giving a bcf file for the reference panel add ALT allele, missing in current format.
    step 1 counts the number of line in bcf header (head variable).
    step 2 from legend file create array with keys line number and values ALT; excluding header (NR>1) from legend file. Then I add the ALT to the bcf but I need to exclude the header of the bcf file. I get the number of lines in "head". Then I start filling at FNR>head and the first ALT is in FNR-head+1. 
    step 3 indexes the bcf file. I am not matching by position and reference as I had entries with same position and reference but different ALT with results in errors. I assume the legend file and bcf are in the same order, which is true, the bcf is transformed from the legend/hap/sample."""
    input:
         lambda wildcards: vcf(config["ref_bcf"])[wildcards.chrom],
         lambda wildcards: vcf(config["ref_legend"])[wildcards.chrom]
    output:
        config['output_dir'] + "/DNA/RP_chr{chrom}_alt_added.bcf" ,
        config['output_dir'] + "/DNA/RP_chr{chrom}_alt_added.bcf.csi"
    shell:
        "head=$(bcftools view -h {input[0]} | wc -l) ; "
        "hm1=$(($head-1)) ; "
        "awk -v head=$head -v hm1=$hm1 "
        " 'FNR==NR{{ if(NR>1) a[NR]=$4;next}}{{if(FNR > head) "
        " $5=a[((FNR - hm1)) ]}}1'  OFS='\t' "
        " <(gzip -dc {input[1]}) <(bcftools view  {input[0]}) "
        " | bcftools view -Ob -o {output[0]}; "
        "tabix {output[0]} "

rule intersect_RP_PEAC:
    """ Extracts and write records from input[0] shared by both input[0] and input[1] using exact allele match. In this case we extract from the reference panel the variants that are present in the PEAC data"""
    input:
        lambda wildcards: config['output_dir'] + "/DNA/RP_chr" + wildcards.chrom + "_alt_added.bcf" ,
        lambda wildcards: config['output_dir'] + "/DNA/PEAC_chr" + wildcards.chrom + "_4PCA.vcf.gz",
        lambda wildcards: config['output_dir'] + "/DNA/RP_chr" + wildcards.chrom + "_alt_added.bcf.csi",
        lambda wildcards: config['output_dir'] + "/DNA/PEAC_chr" + wildcards.chrom + "_4PCA.vcf.gz.tbi"
    output:
        config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz",
        config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz.tbi"
    shell:
         "bcftools isec -n=2 -w1 {input[0]} {input[1]} -Oz -o {output[0]} ; "
         "tabix {output[0]} "

rule intersect_PEAC_RP:
    """ Extracts and write records from input[0] shared by both input[0] and input[1] using exact allele match. Same as above but in reverse order, I just want to make suere files are compatible even if PEAC was imputed with this reference panel, though I dont know if the files used were the same. """
    input:
        lambda wildcards: config['output_dir'] + "/DNA/PEAC_chr" + wildcards.chrom + "_4PCA.vcf.gz",
        lambda wildcards: config['output_dir'] + "/DNA/RP_chr" + wildcards.chrom + "_sub.vcf.gz" ,   
        lambda wildcards: config['output_dir'] + "/DNA/PEAC_chr" + wildcards.chrom + "_4PCA.vcf.gz.tbi",
        lambda wildcards: config['output_dir'] + "/DNA/RP_chr"+ wildcards.chrom + "_sub.vcf.gz.tbi"
    output:
        config['output_dir'] + "/DNA/PEAC_chr{chrom}_sub.vcf.gz"
    shell:
         "bcftools isec -n=2 -w1 {input[0]} {input[1]} -Oz -o {output} "


rule vcf_gds:
    """ From vcf input per chromosome convert to gds while merging into 1 file. I apply it to RP and PEAC""" 
    input:
        expand(config['output_dir'] + "/DNA/PEAC_chr{chrom}_sub.vcf.gz", chrom=vcf(config["geno_vcf"]).keys()) ,
        expand(config['output_dir'] + "/DNA/RP_chr{chrom}_sub.vcf.gz", chrom=vcf(config["geno_vcf"]).keys())
    params:
        method="biallelic.only",
        ##prefix=config['output_dir'] + "/DNA/"
    output:
        peac=config['output_dir'] + "/DNA/PEAC_PCA.gds",
        rp=config['output_dir'] + "/DNA/RP_PCA.gds"
    threads: 2
    script:
        "Rscripts/vcf_PCA2gds.R"

rule RP_PCA:
    """ Compute PCA for reference panel, the matrix of loadings to apply to PEAC samples and applies the matrix of loadings to samples"""
    input:        
        RP=config['output_dir'] + "/DNA/RP_PCA.gds"
    output:
        PC=config['output_dir'] + "/DNA/RP_pcs.rds",
        Loads=config['output_dir'] + "/DNA/RP_loads.rds"
    params:
        ld=0.01,
        maf=0.05,
        method="corr"
    threads: 16
    script:
        "Rscripts/PCA.R"

rule PCs_PEER:
    """ Compute PCs from peac samples by projecting into the reference panel and computes PEER factors based on cqn with GC correction. Prepare inputs for MatQTL"""
    input:
        expr=config['output_dir'] + "/RNA_counts/groups/Genentech_Paired_Baseline_Synovium_RA_Geno.txt",
        peacdata=config['sample_meta'],
        peacGds=config['output_dir'] + '/DNA/PEAC_PCA.gds',
        Loads=config['output_dir'] + "/DNA/RP_loads.rds",
        gene_coord=config['output_dir'] + "/gene_inputs/gene_coord.txt"
    params:
        ld=0.5,
        maf=0.05,
        Nfactors=config['N factors'],
        prefix=["pcs", "peerCqn"]
    output:
        covars=expand(config['output_dir'] + "/matqtl/inputs/{pcs}.{peer}.txt",
                      pcs=["pcs" + str(x) for x in range(1,int(config['N factors'])+1)],
            peer=["peerCqn" + str(x) for x in range(1,int(config['N factors'])+1)]),
        covfix=expand(config['output_dir'] + "/matqtl/inputs/{pcs}.covSexBatch.txt",
                     pcs=["pcs" + str(x) for x in range(1,int(config['N factors'])+1)]),        
        geno=config['output_dir'] + "/matqtl/inputs/genotype.txt",
        expressionCqn=config['output_dir'] + "/matqtl/inputs/gene_expression_cqn.txt",
        geneLoc=config['output_dir'] + "/matqtl/inputs/gene_location.txt",
        snpLoc=config['output_dir'] + "/matqtl/inputs/snp_location.txt"        
    script:
        "Rscripts/PCsPEER.R"

        
        
rule Matqtl:
    """ Compute cis-QTL analysis by linear model. Uses 1:10 PCS and PEER factors"""
    input:
        covars=config['output_dir'] + "/matqtl/inputs/{pcs}.{peerBSex}.txt",     
        geno=config['output_dir'] + "/matqtl/inputs/genotype.txt",
        expressionCqn=config['output_dir'] + "/matqtl/inputs/gene_expression_cqn.txt",
        geneLoc=config['output_dir'] + "/matqtl/inputs/gene_location.txt",
        snpLoc=config['output_dir'] + "/matqtl/inputs/snp_location.txt"             
    params:
        pvOutputThreshold_cis=0.02,
        cisDis=5e5
    output:
        results=config['output_dir'] + "/matqtl/output/{pcs}.{peerBSex}.txt"
    script:
        "Rscripts/matqtl.R"
    
rule Matqtl2:
    """ Compute cis-QTL analysis by linear model. Uses no PCs, within R
    code removes the first row of covariates file which corresponds to the first PC"""
    input:
        covars=config['output_dir'] + "/matqtl/inputs/pcs1.{peerBSex}.txt",     
        geno=config['output_dir'] + "/matqtl/inputs/genotype.txt",
        expressionCqn=config['output_dir'] + "/matqtl/inputs/gene_expression_cqn.txt",
        geneLoc=config['output_dir'] + "/matqtl/inputs/gene_location.txt",
        snpLoc=config['output_dir'] + "/matqtl/inputs/snp_location.txt"             
    params:
        pvOutputThreshold_cis=0.02,
        cisDis=5e5
    output:
        results=config['output_dir'] + "/matqtl/output/pcs0.{peerBSex}.txt"
    script:
        "Rscripts/matqtl2.R"
    
        




        

#PEAC=config['output_dir'] + "/DNA/PEAC_PCA.gds",        
        
# for x in wildcards.group:  ### need to make this rule by gene instead of by chrom ##### use dic{chr:gene}
##### Need to link HospitalNumber (inidivual id) to QMUL id (vcf sample name) #####
#     if x.endswith("Geno"):
#         rule neg_binom_gt:
#             """ Run negative binomial model for samples with genotype grouped by Batch (Ambry or Genentech), Tissue, Timepoint, Reads??, Diagnosis"""
#             input:
#                 lambda wildcards: config['output_dir'] + "/RNA_counts/groups/" + wildcards.group[x] + ".txt" ,
#                 lambda wildcards: config['output_dir'] + "/RNA_counts/groups/" + wildcards.group[x] + "_lib_size.rds" ,
#                 lambda wildcards: vcf(config["geno_vcf"])[wildcards.chrom]
#             output:
#                 config['output_dir'] + "/eqtl/neg_binom_gt/{group}_chr{chrom}.txt"
#             script:
#                 "Rscripts/neg_binom_gt.R"


## make dic{chr: gene list}

## reduce(ebg)
        

## snakemake -j 100 -k --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} --output {cluster.error} -J {cluster.job} "

#snakemake -p -n q # testing (printing) commands, can add --quiet
#snakemake -R somerule ## re run from somerule downstream
