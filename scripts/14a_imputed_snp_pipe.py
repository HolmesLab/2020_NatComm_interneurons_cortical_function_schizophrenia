#!/bin/python

import os
import glob
import numpy as np
import pandas as pd
import subprocess
os.chdir('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn')
from scripts_v2.utilities import writeSlurm, submitSlurm


# Set up paths/directories
proj_dir    = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
out_dir     = os.path.join(proj_dir, 'data/ukb/hsq')
geno_dir    = os.path.join(proj_dir, 'data/ukb/genotyped')
imputed_dir = '/gpfs/milgram/data/UKB/ukb_snp'


# ldsc directories
ldsc_ref_dir = os.path.join(proj_dir, 'data/ukb/imputed/ldsc_ref')
ldsc_dat_dir = os.path.join(proj_dir, 'data/ukb/imputed/ldsc_part')


# external software
plink    = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_1.09/plink'
plink2   = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_2.00/plink2'
ldsc_dir = os.path.join(proj_dir, 'external/ldsc')
ldsc_dir = os.path.join(proj_dir, 'gcta_1.93.0beta/gcta64')


# Step 1

# identify imputed SNPs with INFO > 0.9
names  = ['var_id','snp_id','bp','a1','a2','maf','r','info']
snp_ct = 0
for chrom in np.arange(1,23):
    print(chrom)
    df = pd.read_table('/gpfs/milgram/data/UKB/ukb_snp/ukb_mfi_chr{}_v3.txt'.format(chrom), names=names)
    df['chr'] = str(chrom)
    df['chr_str'] = df['chr'].map(str)

    df['chr_str'] = df['chr_str'].apply(lambda x: x.zfill(2))
    df['snp_id'] = df['chr_str'].map(str) + ':' + df['bp'].map(str)

    # remove duplicated
    df = df.loc[df.duplicated('snp_id') == False]

    df_filt = df.loc[df.a1.isin(['G', 'C', 'T', 'A'])]
    df_filt = df_filt.loc[df_filt.a2.isin(['G', 'C', 'T', 'A'])]
    df_info = df_filt.loc[df_filt['info'] > 0.9]

    snp_ct+=df_info.shape[0]

    # snp positions
    df_info[['chr_str', 'bp', 'bp']].to_csv(os.path.join(proj_dir, 'data/ukb/imputed/snp_lists/ukb_mfi_chr{}_info90_plink_positions.txt'.format(chrom)), sep='\t', index=None, header=None)




# Step 2

# Preprocess imputed UKB bgen data
maf=0.01
hwe=1e-6
mind=0.1
geno=0.1

# partition slurm dir
slurm_out_dir = os.path.join(proj_dir, 'data/ukb/imputed/slurm')
sub_keep_file = os.path.join(proj_dir, 'data/ukb/imputed/sublist_n9713_plink.txt')

# make qctool command and submit to slurm
for chr in np.arange(1,23):
    print(chr)

    # raw ukb bgen/sample paths
    cur_bgen    = '/gpfs/milgram/data/UKB/ukb_snp/ukb_imp_chr' + str(chr) + '_v3.bgen'
    cur_sample  = '/gpfs/milgram/data/UKB/ukb_snp/ukb25163_imp_chr' + str(chr) + '_v3_s487324.sample'

    # create path for the preprocessed chr data
    bgen_out    = os.path.join(proj_dir, 'data/ukb/imputed/snp_data/ukb_imp_chr{}_v3_info090_maf{}_hwe{}_mind{}_geno{}'.format(chr, maf, hwe, mind, geno))
    bgen_out    = bgen_out.replace('0.0', '0').replace('.', '')

    # previously generated SNP file with high INFO variants
    snp_file    = os.path.join(proj_dir, 'data/ukb/imputed/snp_lists/ukb_mfi_chr{}_info90_plink_positions.txt'.format(chr))

    # Command 1: QC bgen files
    plink_cmd   = (f'''{plink2} --bgen {cur_bgen} --sample {cur_sample} --max-alleles 2 --maf {maf} --hwe {hwe} --mind {mind} --geno {geno} -export bgen-1.2 'bits=8' --out {bgen_out} --keep {sub_keep_file} --extract 'range' {snp_file}''')

    # Command 2: convert bgen to bed/bim/fam
    plink_cmd_2 = (f'''{plink2} --bgen {bgen_out}.bgen --sample {bgen_out}.sample --make-bed --out {bgen_out}''')

    # submit both at once
    plink_submit = plink_cmd + '\n\n' + plink_cmd_2

    slurm_file  = os.path.join(slurm_out_dir, 'plink_chr{}'.format(chr))
    cmd_path    = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=18, cmd=plink_submit, stime='06:00:00', jobName=str(chr))
    job_id      = submitSlurm(cmd_path, dependencies=None)




# Step 3: fastGWA

# write chr list for GWAS
mbfile      = os.path.join(proj_dir, 'data/ukb/imputed/chr_files_for_gcta_linear_fastGWA.txt')
mbfile_list = [os.path.join(proj_dir, 'data/ukb/imputed/snp_data/ukb_imp_chr{}_v3_info090_maf01_hwe1e-06_mind01_geno01'.format(i)) for i in np.arange(1,23)]
pd.Series(mbfile_list).to_csv(mbfile, index=None)

# output directory
gwas_out = os.path.join(proj_dir, 'data/ukb/imputed/rsfa_gwas')


# files required to run GCTA fastGWA
pheno  = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_400parcel.txt')
qcovar = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA/RSFA_maf05_brit_n9713_qcovar.txt')
ccovar = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_ccovar.txt')

nthreads = 24
for parcel in np.arange(1,401):
    out_string = 'rsfa_brit_n9713_parcel_{}_400'.format(str(parcel).zfill(3))
    out_file   = os.path.join(gwas_out, out_string)
    gcta_cmd   = f'''{gcta} --mbfile {mbfile} --fastGWA-lr --maf 0.01 --pheno {pheno} --mpheno {parcel} --qcovar {qcovar} --covar {ccovar} --threads {nthreads} --out {out_file}'''
    print(gcta_cmd)

    slurm_file = os.path.join(slurm_out_dir, out_string)
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=nthreads, cmd=gcta_cmd, stime='06:00:00',
                          jobName='gwas_p{}'.format(parcel))
    job_id = submitSlurm(cmd_path, dependencies=None)


# fastGWA GWAS for cluster RSFA estimates
pheno  = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_cluster7_noscale.txt')
for parcel in np.arange(1,8):
    out_string = 'rsfa_brit_n9713_cluster_{}_7'.format(str(parcel).zfill(3))
    out_file   = os.path.join(gwas_out, out_string)
    gcta_cmd   = f'''{gcta} --mbfile {mbfile} --fastGWA-lr --maf 0.01 --pheno {pheno} --mpheno {parcel} --qcovar {qcovar} --covar {ccovar} --threads {nthreads} --out {out_file}'''
    print(gcta_cmd)

    slurm_file = os.path.join(slurm_out_dir, out_string)
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=nthreads, cmd=gcta_cmd, stime='06:00:00',
                          jobName='gwas_p{}'.format(parcel))
    job_id = submitSlurm(cmd_path, dependencies=None)


# Step 4:
#
# Munge Sumstats (i.e. prepare for LDSC)
gwas_out        = os.path.join(proj_dir, 'data/ukb/imputed/rsfa_gwas')
rsfa_gwas_files = glob.glob(os.path.join(gwas_out, 'rsfa_brit_n9713_parcel_*_400.fastGWA'))
rsfa_clust_gwas_files = glob.glob(os.path.join(gwas_out, 'rsfa_brit_n9713_cluster_*_7.fastGWA'))
rsfa_gwas_files = rsfa_gwas_files + rsfa_clust_gwas_files
rsfa_gwas_files.sort()


# LDSC precomputed data is based on hapmap3 snps
hapmap_snps = os.path.join(ldsc_ref_dir, 'w_hm3.snplist')

munge_cmd = os.path.join(ldsc_dir, 'munge_sumstats.py')
for gwas in rsfa_gwas_files:
    print(gwas)
    
    # make LDSC munge command
    ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/mdd_sst/external/ldsc/ \n'
    ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
    ldsc_cmd = ldsc_cmd + f'''{munge_cmd} --sumstats {gwas} --merge-alleles {hapmap_snps} --out {gwas}'''
    
    # submit to cluster
    slurm_file = os.path.join(slurm_out_dir, 'munge_'+gwas.split('/')[-1] )
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=2, cmd=ldsc_cmd, stime='06:00:00',
                          jobName=gwas.split('/')[-1])
    job_id = submitSlurm(cmd_path, dependencies=None)



# Step 5:
#
# LDSC heritability
ld_ref  = 'eur_w_ld_chr'
ldsc_py = os.path.join(ldsc_dir, 'ldsc.py')
for gwas in rsfa_gwas_files:
    munge_gwas = gwas + '.sumstats.gz'
    gwas_out   = munge_gwas.replace('.sumstats.gz','_' + ld_ref + '_ldsc_h2').replace('.','_')
    print(gwas_out)

    # reference LD info 1000G
    ref_w_ld = os.path.join(ldsc_ref_dir, ld_ref + '/')
    
    # build LDSC heritability command
    ldsc_cmd = 'cd /gpfs/milgram/project/holmes/kma52/mdd_sst/external/ldsc/ \n'
    ldsc_cmd = ldsc_cmd + 'source activate ldsc \n\n'
    ldsc_cmd = ldsc_cmd + f'''{ldsc_py} --h2 {munge_gwas} --ref-ld-chr {ref_w_ld} --w-ld-chr {ref_w_ld} --out {gwas_out}'''
    print(ldsc_cmd)

    slurm_file = os.path.join(slurm_out_dir, 'ldsc_h2_'+gwas.split('/')[-1] )
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=2, cmd=ldsc_cmd, stime='06:00:00',
                          jobName=gwas.split('/')[-1])
    job_id = submitSlurm(cmd_path, dependencies=None)








# end


