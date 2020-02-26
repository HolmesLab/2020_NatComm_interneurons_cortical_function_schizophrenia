#!/bin/python

import os
import pandas as pd
import subprocess
os.chdir('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn')
from scripts_v2.utilities import writeSlurm, submitSlurm


# external software
plink = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_1.09/plink'
gcta  = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/gcta_1.91.1beta/gcta64'


# Set up paths/directories
proj_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
out_dir  = os.path.join(proj_dir, 'data/ukb/hsq')
geno_dir = os.path.join(proj_dir, 'data/ukb/genotyped')


# func phenotype (RSFA)
n = '9713'
analysis_dir = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA')


# make slurm dir
analysis_dir = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA')
slurm_dir    = os.path.join(analysis_dir, 'slurm')
if not os.path.exists(slurm_dir):
    os.mkdir(slurm_dir)

# hsq dir
hsq_dir = os.path.join(analysis_dir, 'hsq')
if not os.path.exists(hsq_dir):
    os.mkdir(hsq_dir)

# hsq partition dir
part_dir = os.path.join(analysis_dir, 'hsq_partition')
if not os.path.exists(part_dir):
    os.mkdir(part_dir)

# for each MAF level
pheno     = os.path.join(analysis_dir, 'RSFA_brit_n'+n+'_cluster7_noscale.txt')
out_dir   = os.path.join(proj_dir, 'data/ukb/hsq')


# GRM
# genotype dat
# side note, the n10313 refers to the number of subjects in the GRM prior to QC and removal of cryptic relatedness
# it actually contains 9,713 individuals
bfile_base = 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025'
bfile      = os.path.join(geno_dir, bfile_base)
gcta_file  = os.path.join(geno_dir, bfile_base)


# GCTA formated RSFA values and genetic/phenotypic covariates
pheno  = os.path.join(analysis_dir, 'RSFA_brit_cluster7_noscale.txt')
qcovar = os.path.join(analysis_dir, 'RSFA_maf05_brit_n'+n+'_qcovar.txt')
ccovar = os.path.join(analysis_dir, 'RSFA_brit_n'+n+'_ccovar.txt')



# overall heritability clusters of each overarching cluster, calculated with GCTA
net_names = ['cingulo_opercular','visual','parietal','motor','prefrontal','limbicA', 'limbicB']
for idx,net in enumerate(net_names):
    print(net)
    gcta_cmd   = f'''{gcta} --grm {gcta_file} --reml --pheno {pheno} --covar {ccovar} --qcovar {qcovar} --mpheno {idx+1} --out {hsq_dir}/{bfile_base}_{net}_clust7_noscale --thread-num 2'''
    slurm_file = os.path.join(slurm_dir, 'RSFA_' + net + '_' + bfile_base + '_heritability')
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=2, cmd=gcta_cmd, stime='06:00:00', jobName=net)
    job1_id    = submitSlurm(cmd_path, dependencies=None)


# Heritability of parcel-wise measures
# schaeffer dir/names
parcel_pheno = os.path.join(analysis_dir, 'RSFA_brit_n'+n+'_400parcel.txt')
parcel_names = pd.read_csv(os.path.join(analysis_dir, 'RSFA_brit_n'+n+'_400parcel_colnames.txt'), header=None)

# GRM
bfile_base = 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025'
bfile      = os.path.join(geno_dir, bfile_base)
gcta_file  = os.path.join(geno_dir, bfile_base)

# overall heritability of RSFA, across each schaefer parcel 
for idx,net in enumerate(parcel_names[0]):
    print(idx)
    print(net)
    gcta_cmd   = f'''{gcta} --grm {gcta_file} --reml --pheno {parcel_pheno} --covar {ccovar} --qcovar {qcovar} --mpheno {idx+1} --out {hsq_dir}/{net}_{bfile_base} --thread-num 2'''
    slurm_file = os.path.join(slurm_dir, 'RSFA_' + net + '_' + bfile_base + '_heritability')
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=2, cmd=gcta_cmd, stime='06:00:00', jobName=net)
    job1_id    = submitSlurm(cmd_path, dependencies=None)



lower = 1
upper = 500

# Partitioned heritability
pvalb_rsid = pd.read_csv(os.path.join(proj_dir,'gene_lists/rsid_lists/rr_ahba_rsid_PVALB_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_'+lower+'_'+upper+'.txt'), header=None)
sst_rsid   = pd.read_csv(os.path.join(proj_dir,'gene_lists/rsid_lists/rr_ahba_rsid_SST_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_'+lower+'_'+upper+'.txt'), header=None)
all_rsids  = pd.concat([pvalb_rsid, sst_rsid])

# create text file containing all SNPs in the first 1-500 sst/pvalb
snp_list = os.path.join(proj_dir,'gene_lists/rsid_lists/rr_SNPlist_PVALBSST_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_'+lower+'_'+upper+'.txt')
all_rsids.to_csv(snp_list, index=False, header=None)

# set up paths
analysis_dir = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA)
slurm_dir    = os.path.join(analysis_dir, 'slurm')
part_dir     = os.path.join(analysis_dir, 'hsq_partition')
bfile_base   = 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025'
gcta_file    = os.path.join(geno_dir, bfile_base)

out_file = os.path.join(part_dir, bfile_base + '_ahba_'+lower+'-'+upper+'_exclude')

# submit to slurm
exclude_cmd = f'''{gcta} --bfile {gcta_file} --make-grm --exclude {snp_list} --out {out_file} --thread-num 8'''
print(exclude_cmd)
slurm_file  = os.path.join(slurm_dir, 'RSFA_ahba_'+lower+'-'+upper+'_exclude')
cmd_path    = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=8, cmd=exclude_cmd, stime='06:00:00', jobName='ahba_'+lower+'-'+upper+'_exclude')
job1_id     = submitSlurm(cmd_path, dependencies=None)

                            
# compute GRMs for the PVALB and SST SNP partitions
for gene in ['PVALB','SST']:
    print(gene)
    snp_list = os.path.join(proj_dir,'gene_lists/rsid_lists/rr_ahba_rsid_'+gene+'_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_'+lower+'_'+upper+'.txt')

    # build output paths
    uniq_id  = gene+'_'+str(lower)+'_'+str(upper)
    out_file = os.path.join(part_dir, bfile_base + '_' + uniq_id + '_extract')
    exclude_out_file = os.path.join(part_dir, bfile_base + '_' + uniq_id + '_exclude')

    # create GCTA commands                        
    extract_cmd = f'''{gcta} --bfile {gcta_file} --make-grm --extract {snp_list} --out {out_file} --thread-num 8 \n\n'''
    extract_cmd = extract_cmd + f'''{gcta} --bfile {gcta_file} --make-grm --exclude {snp_list} --out {exclude_out_file} --thread-num 8'''
    print(extract_cmd)
                            
    # submit to slurm
    slurm_file  = os.path.join(slurm_dir, 'RSFA_'+gene+'__'+lower+'_'+upper+'_make_grm')
    cmd_path    = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=8, cmd=extract_cmd, stime='06:00:00', jobName='extract')
    job1_id     = submitSlurm(cmd_path, dependencies=None)



# make the GRM list (PVALB/SST/ALL)
grm_list = list()
grm_list.append(os.path.join(part_dir, 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_PVALB_'+lower+'_'+upper+'_extract'))
grm_list.append(os.path.join(part_dir, 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_SST_'+lower+'_'+upper+'_extract'))
grm_list.append(os.path.join(part_dir, 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_ahba_'+lower+'-'+upper+'_exclude'))
grm_list_path = os.path.join(part_dir, 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_grm_list_'+lower+'-'+upper+'.txt')
pd.DataFrame(grm_list).to_csv(grm_list_path, index=False, header=None)


# cluster partitioned heritability clusters
analysis_dir = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA')
part_dir     = os.path.join(analysis_dir, 'hsq_partition')
slurm_dir    = os.path.join(analysis_dir, 'slurm')
                            
# paths to phenotype and covariate information
pheno  = os.path.join(analysis_dir, 'RSFA_brit_n9713_cluster7_noscale.txt')
qcovar = os.path.join(analysis_dir, 'RSFA_maf05_brit_n9713_qcovar.txt')
ccovar = os.path.join(analysis_dir, 'RSFA_brit_n9713_ccovar.txt')

# GCTA based partitioned heritability for each of the 7 clustered RSFA values 
net_names = ['cingulo-opercular', 'visual', 'parietal', 'motor', 'limbicB', 'prefrontal','limbicA']
for idx, net in enumerate(net_names):
    print(idx)
    gcta_cmd   = f'''{gcta} --mgrm {grm_list_path} --reml --pheno {pheno} --covar {ccovar} --qcovar {qcovar} --mpheno {idx + 1} --out {part_dir}/{bfile_base}_{net}_partitioned_clust7_{lower}-{upper} --thread-num 4'''
    slurm_file = os.path.join(slurm_dir, 'RSFA_' + net + '_' + bfile_base + '_heritability_'+lower+'-'+upper+'')
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=4, cmd=gcta_cmd, stime='06:00:00', jobName=net)
    job1_id    = submitSlurm(cmd_path, dependencies=None)



                            
# partitioned heritability of RSFA in each of the 400 Schaefer parcels 
analysis_dir  = os.path.join(proj_dir, 'data/ukb/rsfa/pheno_RSFA')
part_dir      = os.path.join(analysis_dir, 'hsq_partition')
slurm_dir     = os.path.join(analysis_dir, 'slurm')
grm_list_path = os.path.join(part_dir, 'ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_grm_list_'+lower+'-'+upper+'.txt')

# paths to phenotypes/covariates 
qcovar       = os.path.join(analysis_dir, 'RSFA_maf05_brit_n9713_qcovar.txt')
ccovar       = os.path.join(analysis_dir, 'RSFA_brit_n9713_ccovar.txt')
parcel_pheno = os.path.join(analysis_dir, 'RSFA_brit_n9713_400parcel.txt')
parcel_names = pd.read_csv(os.path.join(analysis_dir, 'RSFA_brit_n9713_400parcel_colnames.txt'), header=None)

# partitioned heritability of individual parcels
for idx, net in enumerate(parcel_names[0]):
    print(idx)
    print(net)
    gcta_cmd   = f'''{gcta} --mgrm {grm_list_path} --reml --pheno {parcel_pheno} --covar {ccovar} --qcovar {qcovar} --mpheno {idx + 1} --out {part_dir}/v2_{net}_{bfile_base}_{lower}_{upper} --thread-num 2'''
    slurm_file = os.path.join(slurm_dir, 'RSFA_' + net + '_' + bfile_base + '_heritability')
    cmd_path   = writeSlurm(slurm_file=slurm_file, partition='short', nthreads=2, cmd=gcta_cmd, stime='06:00:00', jobName=net)
    job1_id    = submitSlurm(cmd_path, dependencies=None)






