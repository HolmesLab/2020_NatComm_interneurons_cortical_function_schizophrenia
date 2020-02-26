#!/bin/bash

# project directory
proj_dir=/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn


# set up paths for MGMA annotation
magma=${proj_dir}/external/magma
scz_gwas_snps=${proj_dir}/data/scz_gwas/snp_loc_ckqny.scz2snpres.txt
genome_build=${proj_dir}/data/scz_gwas/NCBI37.3.gene.loc
scz_out=${proj_dir}/data/scz_gwas/snp_loc_ckqny.scz2snpres


# annotate genes
${magma} --annotate window=5,5 \
--snp-loc ${scz_gwas_snps} \
--gene-loc ${genome_build} \
--out ${scz_out}


# Run MAGMA
ref_bfile=/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/scz_gwas/g1000_eur
scz_gwas=/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/scz_gwas/ckqny.scz2snpres

# gene annotation for the Ripke SCZ GWAS
${magma}  \
--bfile ${ref_bfile} \
--pval ${scz_gwas} N=150064 \
--gene-annot ${scz_out}.genes.annot \
--out ${scz_out}




