#!/bin/python

import os
import glob
import numpy as np
import pandas as pd
import subprocess
os.chdir('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn')
from scripts_v2.utilities import writeSlurm, submitSlurm


# set up directories
base_dir     = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ldsc_dat_dir = os.path.join(base_dir, 'data/ukb/imputed/ldsc_part')
ldsc_ref_dir = os.path.join(base_dir, 'data/ukb/imputed/ldsc_ref')


# paths to previously generated gene location references file (generated from biomaRt)
gene_coord_file = os.path.join(ldsc_dat_dir, 'ahba_hg37_gene_coords.txt')


# Step 1: create annotation files for SST/PVALB top 500 gene sets
# ---------
for gene in ['SST','PVALB']:

    # create LDSC annotation file for the top 500 SST/PVALB genes
    annot_base = os.path.join(ldsc_dat_dir, 'GCTA_{}_ahba_top500_ensembl'.format(gene))
    for chr in np.arange(1,23):

        # required LDSC data paths
        annot_file      = '{}_chr{}'.format(annot_base, chr)
        hapmap_file     = os.path.join(ldsc_ref_dir, 'hapmap3_snps/hm.{}.snp'.format(chr))
        eur_1000G_file  = os.path.join(ldsc_ref_dir, '1000G_Phase1_plinkfiles/1000G_plinkfiles/1000G.mac5eur.{}'.format(chr))

        # commands for LDSC annotation file generation
        ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/ldsc/ \n'
        ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
        ldsc_cmd = ldsc_cmd + f'''python make_annot.py --gene-set-file {annot_base}.txt --gene-coord-file {gene_coord_file} --windowsize 10000 --bimfile {eur_1000G_file}.bim --annot-file {annot_file}.annot.gz \n\n'''
        ldsc_cmd = ldsc_cmd + f'''python ldsc.py --l2 --bfile {eur_1000G_file} --ld-wind-cm 1 --annot {annot_file}.annot.gz --thin-annot --out {annot_file} --print-snps {hapmap_file} \n\n'''
        print(ldsc_cmd)

        # submit to cluster
        slurm_file = os.path.join(ldsc_dat_dir, 'pvalb_order_ldsc_files/GCTA_{}_ahba_top500_ensembl_chr{}'.format(gene, chr))
        cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=3, cmd=ldsc_cmd, stime='06:00:00', jobName=str(chr))
        submitSlurm(cmd_path, dependencies=None)



# Step 2: partitioned heritability
# ---------
gwas_out = os.path.join(base_dir, 'data/ukb/imputed/rsfa_gwas')
rsfa_clust_gwas_files = glob.glob(os.path.join(gwas_out, 'rsfa_brit_n9713_cluster_*_7.fastGWA.sumstats.gz'))
rsfa_clust_gwas_files.sort()

# ldsc ref files
baseline_files = os.path.join(ldsc_ref_dir, '1000G_Phase1_baseline_ldscores/baseline/baseline.')
frq_files      = os.path.join(ldsc_ref_dir, '1000G_frq/1000G.mac5eur.')
weight_files   = os.path.join(ldsc_ref_dir, 'weights_hm3_no_hla/weights.')

for gene in ['SST','PVALB']:
    annot_base = os.path.join(ldsc_dat_dir, 'GCTA_{}_ahba_top500_ensembl'.format(gene))
    annot_ldsc = os.path.join(annot_base + '_chr')

    for rsfa_gwas in rsfa_clust_gwas_files:
        print(rsfa_gwas)
        rsfa_name = rsfa_gwas.split('/')[-1]
        ldsc_out  = annot_base + '_' + rsfa_name.split('.')[0]

        # commands for LDSC annotation file generation
        ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/ldsc/ \n'
        ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
        ldsc_cmd = ldsc_cmd + f'''python ldsc.py --h2 {rsfa_gwas} --w-ld-chr {weight_files} --ref-ld-chr {annot_ldsc},{baseline_files} --overlap-annot --frqfile-chr {frq_files} --out {ldsc_out} --print-coefficients\n\n'''
        print(ldsc_cmd)

        cmd_path   = writeSlurm(slurm_file=ldsc_out, partition='short', nthreads=3, cmd=ldsc_cmd, stime='06:00:00', jobName=str(chr))
        submitSlurm(cmd_path, dependencies=None)



# do the same partitioned heritability analysis, but individually for each RSFA parcel
rsfa_gwas_files = glob.glob(os.path.join(gwas_out, 'rsfa_brit_n9713_parcel_*_400.fastGWA.sumstats.gz'))
rsfa_gwas_files.sort()

for gene in ['SST','PVALB']:
    annot_base = os.path.join(ldsc_dat_dir, 'GCTA_{}_ahba_top500_ensembl'.format(gene))
    annot_ldsc = os.path.join(annot_base + '_chr')

    for rsfa_gwas in rsfa_gwas_files:
        print(rsfa_gwas)

        rsfa_name = rsfa_gwas.split('/')[-1]
        ldsc_out  = annot_base + '_' + rsfa_name.split('.')[0]

        # commands for LDSC annotation file generation
        ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/ldsc/ \n'
        ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
        ldsc_cmd = ldsc_cmd + f'''python ldsc.py --h2 {rsfa_gwas} --w-ld-chr {weight_files} --ref-ld-chr {annot_ldsc},{baseline_files} --overlap-annot --frqfile-chr {frq_files} --out {ldsc_out} --print-coefficients\n\n'''
        print(ldsc_cmd)

        cmd_path   = writeSlurm(slurm_file=ldsc_out, partition='short', nthreads=3, cmd=ldsc_cmd, stime='06:00:00', jobName=str(chr))
        submitSlurm(cmd_path, dependencies=None)



