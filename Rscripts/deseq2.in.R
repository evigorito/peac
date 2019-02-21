source("/home/ev250/emedlab/snake-pipe/Rfunctions/inputs.eQTL.dseq2.R")

in.deseq2(gene=snakemake@wildcards[['gene']],
          chr=as.numeric(snakemake@params[['chrom']]),
          snps=as.numeric(snakemake@params[['snps']]),
          counts.f=snakemake@input[['counts']],
          gene.coord=snakemake@input[['genecoord']],
          vcf=snakemake@input[['vcf']],
          nhets=as.numeric(snakemake@params[['nhets']]),
          tag.threshold=as.numeric(snakemake@params[['tag']]),
          out=snakemake@params[['out']],
          missing=as.numeric(snakemake@params[['missing']])
          )
