## Snakefile for emedlab


## https://snakemake.readthedocs.io/en/stable/project_info/faq.html#i-want-to-configure-the-behavior-of-my-shell-for-all-rules-how-can-that-be-achieved-with-snakemake

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
     <unique_sample_id> <unique_individual_id><Batch><Reads><Timepoint><Tissue><Diagnosis><GenotypeDir><fastq1.R1_path> <fastq1.R2_path> <fastq2.R1_path> <fastq2.R2_path>. This function produces a dictionary of sample_id keys and values(individualID, Batch, Reads, Timepoint, Tissue, Diagnosis, Dir with genotype, and a sub-list of (fq1, fq2, ...,) 
     Input file prepared in /home/ev250/emedlab/snake-pipe/Rscripts/files4fastqc.R
    """
    with open(config['sample_file'], "r") as f:
        samp_dict = {}
        for line in f:
            words = line.strip().split(",")
            words2 = words[8:]
            while 'NA'in words2:
                words2.remove('NA')
            words = words[:8]
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
    sub_values = [v[1:7] for k, v in samp_dict.items()]

    ## get unique combinations of variables
    u_val = [list(x) for x in set(tuple(x) for x in sub_values)]
    ## Prepare dic with key: Batch_reads_tiempoint_tissue_diagnosis_geno, batch Ambry or genentech. Only one sample used if more of the same type available for the same indidual
    genen_u_ind = group_helper(samp_dict,u_val,"Genentech")
    ambry_u_ind = group_helper(samp_dict,u_val,"Ambry")
    return(genen_u_ind.update(ambry_u_ind))


def group_helper(samp_dict, u_val, batch):
    """ Helper function to iterate over Genentech or Ambry. Creates a dictionary with keys Batch_reads_tiempoint_tissue_diagnosis_geno, batch Ambry or genentech, and values sample_name as in fasq
 """
    u_val = [list(x) for x in set(tuple(x[1:]) for x in u_val if x[0].startswith(batch))]
    u_ind = {}
    for x in range(0,len(u_val)):
    ## prepare values
        val = [list(v) for v in set(tuple(v[:7]) for (k,v) in samp_dict.items() if v[2:7]==u_val[x] if v[1].startswith(batch))]
        val = [k for (k,v) in samp_dict.items() for x in val if v[:7] == x ]
        ## prepare key
        key= u_val[x].copy()
        if key[4] != "NA":
            key[4] = "Geno"
        else:
            key.remove("NA")
            
        key.insert(0,batch)
        key = "_".join(key)
        u_ind[key] = val
    return(u_ind)


SAMPLE = read_samples().keys()
#READ = [item[2] for item in read_samples().values()]
#Read = set(READ)

rule all:
    input:
        ##expand(config['output_dir'] + "/STAR/{read}/{sample}/Aligned.sortedByCoord.out.bam" , zip, sample=read_samples().keys(), read=[item[2] for item in read_samples().values()])
        expand(config['output_dir'] + "/RNA_counts/{sample}.txt" , sample=SAMPLE)
        expand(config['output_dir'] + "/RNA_counts/{group}.txt" , group=group_samples.keys() )

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
    """ Get exons per gene from the gft annotation file used for the alignment, it will be used to get total raw counts per gene """
    input:
        config['ref_gtf']
    output:
        config['ebg']
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
    """ Group total gene counts by tissue, timepoint, diagnosis and batch"""
    input:
        expand(config['output_dir'] + "/RNA_counts/{sample}.txt" , sample=SAMPLE)
    output:       
        config['output_dir'] + "/RNA_counts/{group}.txt"
    script:
         "Rscripts/group_gene_counts.R"





         
    # output:
    #     config['ebg'] + "/" + config['built'] + {Read}_ebg.rds
        
        
# ## snakemake --use-conda --cores 1
## snakemake -j 20 -k --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.N} -n {cluster.n}  -t {cluster.time} --output {cluster.error} -J {cluster.job} "
