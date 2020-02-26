#!/bin/bash

# set up paths to required software
source /gpfs/milgram/project/holmes/HOLMES_UKB/projects/HOLMES_UKB_init.bash

# set up directories
proj_dir=/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn
out_dir=${proj_dir}/data/ukb/genotyped
geno_dir=/gpfs/milgram/data/UKB/REPOSITORY/GWAS/NON_IMPUTED


# plink/gcta paths
plink2=/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_2.00/plink2
plink=/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_1.09/plink
gcta=/gpfs/milgram/project/holmes/HOLMES_UKB/external/gcta_1.91.1beta/gcta64


# subjects with rsfc data
n='10313'
sub_list=${proj_dir}/data/ukb/genotyped/ukb_imaging_subs_before_cryptic_relatedness_n${n}.txt
num_subs=`cat ${sub_list} | wc -l`


# perform SNP censoring seperately for each CHR
for maf in '01' '05';
do
    echo $maf
    for chr in {1..22};
    do
        snp_out=${out_dir}/ukb_geno_chr${chr}_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.${maf}_n${n}
        slurm_file=${out_dir}/slurm/chr${chr}_maf${maf}_preprocess_snp.txt
        slurm_out=${out_dir}/slurm/chr${chr}_maf${maf}_preprocess_snpOut.txt
        echo ${chr}
        cat > ${slurm_file} << EOF
#!/bin/bash
#SBATCH --partition=short
#SBATCH --output=${slurm_out}
#SBATCH --job-name=chr${chr}
#SBATCH --ntasks=1 --cpus-per-task=4 --nodes=1
#SBATCH --time=6:00:00
source /gpfs/milgram/project/holmes/HOLMES_UKB/projects/HOLMES_UKB_init.bash
${plink2} --bfile ${geno_dir}/ukb_cal_chr${chr}_v2 \
                      --keep ${sub_list} \
                      --make-bed \
                      --geno 0.02 \
                      --mind 0.1 \
                      --hwe 1e-6 \
                      --maf 0.${maf} \
                      --out ${snp_out}
EOF
        chmod 770 ${slurm_file}
        sbatch < ${slurm_file}
    done
done


# create file list to merge genotypes
for maf in '01' '05';
    do
    merge_file=${out_dir}/ukb_imp_v3_caucasian_merge_list_n${n}.txt
    # initiate file with chr 1
    echo ${out_dir}/ukb_geno_chr1_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.${maf}_n${n} > ${merge_file}
    # add chr 2-22
    for chr in {2..22}; do
      snp_out=${out_dir}/ukb_geno_chr${chr}_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.${maf}_n${n}
      echo ${snp_out} >> ${merge_file}
    done

    # combine individual chromosomes
    bfile_comb=${out_dir}/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.${maf}_n${n}
    ${plink} --merge-list ${merge_file} \
             --make-bed \
             --out ${bfile_comb}
done


# combine data file names
bfile_comb_01=${out_dir}/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.01_n${n}
bfile_comb_05=${out_dir}/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n${n}


# create GRM from the SNP data
${gcta} --bfile ${bfile_comb_05} --autosome --make-grm --out ${bfile_comb_05} --thread-num 4


# threshold GRM
${gcta} --grm ${bfile_comb_05} --grm-cutoff 0.025 --make-grm --out ${bfile_comb_05}_rm025

for chr in {1..22};
do
    orig_bfile=${out_dir}/ukb_geno_chr${chr}_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n${n}
    new_bfile=${out_dir}/ukb_geno_chr${chr}_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n${n}_rm025
    ${plink2} --bfile ${orig_bfile} \
                          --keep ${bfile_comb_05}_rm025.grm.id \
                          --make-bed \
                          --out ${new_bfile}
done

# remove cryptic relatedness
${plink2} --bfile ${bfile_comb_05} \
                      --keep ${bfile_comb_05}_rm025.grm.id \
                      --make-bed \
                      --out ${bfile_comb_05}_rm025


# use the unrelated subjects from MAF=05 in order to match samples
# remove cryptic relatedness
${plink2} --bfile ${bfile_comb_01} \
                      --keep ${bfile_comb_05}_rm025.grm.id \
                      --make-bed \
                      --out ${bfile_comb_01}_rm025

${gcta} --bfile ${bfile_comb_01}_rm025 --autosome --make-grm --out ${bfile_comb_01}_rm025 --thread-num 4


# Calculate top 10 PCA eigenvectors
${plink} --bfile ${bfile_comb_05}_rm025 \
     --pca 10    \
     --out ${bfile_comb_05}_rm025


${plink} --bfile ${bfile_comb_01}_rm025 \
     --pca 10    \
     --out ${bfile_comb_01}_rm025


#end
