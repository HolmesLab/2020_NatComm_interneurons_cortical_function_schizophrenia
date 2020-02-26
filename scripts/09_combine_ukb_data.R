library(tidyverse)
library(dplyr)


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# read UKB phenotype data
ukb_df = read_csv(paste0(base_dir, '/data/ukb/ukb_pheno_processed.csv'))


# read raw RSFA data
rsfa_df  = read_csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_MNI_filtered_func_data_clean_despike_blur4_noGS_noWM_noCSF_schaeff400.csv'))
colnames(rsfa_df) = gsub('7Networks', 'RSFA_Net7', colnames(rsfa_df))


# combine imaging df
imaging_data = rsfa_df

# remove subjects with no imaging data
rsfa_dat = imaging_data %>% select(-one_of('UKB_ID')) # remove ukb id column


# filter subjects with missing data
imaging_data = imaging_data[which(rowSums(rsfa_dat == 0) == 0),]


# Merge pheno/RSFA data and remove subjects with motion greater than 0.20
ukb_imaging = merge(x=ukb_df, y=imaging_data, by.x='f.eid', by.y='UKB_ID')


# write data
write_csv(x=ukb_imaging, path=paste0(base_dir, '/data/ukb/ukb_imaging_beforeQC.csv'))
ukb_imaging = read_csv(paste0(base_dir, '/data/ukb/ukb_imaging_beforeQC.csv'))


ukb_imaging = ukb_imaging %>% filter(MRI_REST_motion.2.0 <= 0.20)

# calculate age at imaging visit
elapsed_months = function(end_date, start_date) {
  ed = as.POSIXlt(end_date)
  sd = as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

# convert month strings to numeric
months = c('January','February','March','April','May','June','July','August','September','October','November','December')
ukb_imaging$monthnum_of_birth.0.0 = match( ukb_imaging$month_of_birth.0.0, months )
ukb_imaging$bday = paste(ukb_imaging$year_of_birth.0.0, ukb_imaging$monthnum_of_birth.0.0, '01', sep='-')

# number of months between DOB and imaging visit/12
image_age = NULL
for (row in 1:nrow(ukb_imaging)){
  age = elapsed_months(ukb_imaging$Date_of_attending_assessment_centre.2.0[row], ukb_imaging$bday[row])/12
  image_age = c(image_age, age)
}
ukb_imaging$image_visit_age = image_age


# blood pressure- calculate average across manual/auto readings, removing NAs where possible
ukb_imaging$diastolic.2.0  = rowMeans(cbind(ukb_imaging$diastolic_auto.2.0, ukb_imaging$diastolic_manual.2.0), na.rm=T)
ukb_imaging$systolic.2.0   = rowMeans(cbind(ukb_imaging$systolic_auto.2.0, ukb_imaging$systolic_manual.2.0), na.rm=T)

# remove sex mismatch -- 279 subjects have NA values for genetic_sex.0.0
ukb_imaging = ukb_imaging[which(ukb_imaging$genetic_sex.0.0 == ukb_imaging$sex.0.0),]

# replace missing BMI values
bmi_handcalc = ukb_imaging$weight_preimaging.2.0/((ukb_imaging$height_standing.0.0/100)^2)
ukb_imaging$BMI.2.0[is.na(ukb_imaging$BMI.2.0)] = bmi_handcalc[is.na(ukb_imaging$BMI.2.0)]
# same for height
ukb_imaging$height_standing = ukb_imaging$height_standing.2.0
ukb_imaging$height_standing[is.na(ukb_imaging$height_standing.2.0)] = ukb_imaging$height_standing.0.0[is.na(ukb_imaging$height_standing.2.0)]
#summary(ukb_df$height_standing)


# variables in the regression (row-wise deletion)
ukb_imaging[['scanner_Y_pos.2.0']] = ukb_imaging[['scanner_Y_pos .2.0']]
regress_vars = c('f.eid','sex.0.0','image_visit_age','BMI.2.0','height_standing','weight_preimaging.2.0',
                  'MRI_REST_motion.2.0','MRI_T1_invSNR.2.0','MRI_REST_invSNR_preproc.2.0', 'diastolic.2.0', 'systolic.2.0',
                  'MRI_T1_vetntCSF_norm.2.0', 'MRI_T1_volGM_WM_headNorm.2.0', 'UK_Biobank_assessment_centre.2.0',
                  'scanner_X_pos.2.0','scanner_Z_pos.2.0','scanner_Y_pos.2.0','scanner_table_pos.2.0')
check_missing = ukb_imaging[regress_vars]
regress_vars[!regress_vars %in% colnames(ukb_imaging)]
# give feedback about which variables have missing data
for (r in regress_vars){
  print(r)
  print(length(which(is.na(ukb_imaging[[r]]))))
}


# white matter lesions
lesion_df   = read_csv('/gpfs/milgram/project/holmes/kma52/ukb_pymood/data/ukb/dataframes/lesion_df.csv')
ukb_imaging = merge(x=ukb_imaging, y=lesion_df, by.x='f.eid', by.y='UKB_ID')


# Subset to only Caucasian subjects
ukb_cauc_gene = ukb_imaging %>% filter(genetic_ethnicity.0.0 == 'Caucasian')


# median value replacement for missing diastolic/systolic BP
ukb_cauc_gene$systolic.2.0[is.na(ukb_cauc_gene$systolic.2.0)] = ukb_cauc_gene$systolic.0.0[is.na(ukb_cauc_gene$systolic.2.0)]
ukb_cauc_gene$systolic.2.0[is.na(ukb_cauc_gene$systolic.2.0)] = median(ukb_cauc_gene$systolic.2.0, na.rm=T)

ukb_cauc_gene$diastolic.2.0[is.na(ukb_cauc_gene$diastolic.2.0)] = ukb_cauc_gene$diastolic.0.0[is.na(ukb_cauc_gene$diastolic.2.0)]
ukb_cauc_gene$diastolic.2.0[is.na(ukb_cauc_gene$diastolic.2.0)] = median(ukb_cauc_gene$diastolic.2.0, na.rm=T)


# remove invSNR outliers
thresh        = mean(ukb_cauc_gene$MRI_REST_invSNR_fullPreproc.2.0) + (4*sd(ukb_cauc_gene$MRI_REST_invSNR_fullPreproc.2.0))
ukb_cauc_gene = ukb_cauc_gene %>% filter(MRI_REST_invSNR_fullPreproc.2.0 <= thresh)


# calc brain size to normalize lesion volume
ukb_cauc_gene$brain_size  = ukb_cauc_gene$MRI_T1_volGM_WM.2.0 + ukb_cauc_gene$MRI_T1_ventCSF.2.0
ukb_cauc_gene$lesion_vol[is.na(ukb_cauc_gene$lesion_vol)] = median(ukb_cauc_gene$lesion_vol, na.rm=T)
ukb_cauc_gene$lesion_norm = as.numeric(scale(resid(lm('lesion_vol ~ brain_size', data=ukb_cauc_gene))))


ggplot(ukb_cauc_gene[ukb_cauc_gene$lesion_norm < 10,], aes(x=lesion_norm))+
  geom_density(color="darkblue", fill="lightblue") +
    geom_vline(xintercept = 2) +
    theme_classic() +
    ggtitle('distribution of brain size normalized lesion volumes') +
    theme(plot.title = element_text(size = 10))

dim(ukb_cauc_gene)


# recode sex to numeric
ukb_cauc_gene$sex_num.0.0 = recode(ukb_cauc_gene$sex.0.0, 'Male'=1, 'Female'=0)
table(ukb_cauc_gene$sex_num.0.0)

# rsfa column names
rsfa_names = colnames(ukb_cauc_gene)[grep('RSFA', colnames(ukb_cauc_gene))]

# calculate average RSFA across 400 parcels
print(dim(ukb_cauc_gene[rsfa_names]))
ukb_cauc_gene$rsfa_mean = rowMeans(ukb_cauc_gene[rsfa_names])

summary(ukb_cauc_gene$rsfa_mean)

ggplot(ukb_cauc_gene, aes(x=rsfa_mean))+
  geom_density(color="darkblue", fill="lightblue") +
    theme_classic() +
    ggtitle('Average RSFA per subject') +
    theme(plot.title = element_text(size = 10))


#
thr = mean(ukb_cauc_gene$rsfa_mean)+4*sd(ukb_cauc_gene$rsfa_mean)
length(which(ukb_cauc_gene$rsfa_mean >= thr))

# remove genetic_heterozygosity outliers
uthr = mean(ukb_cauc_gene$genetic_heterozygosity.0.0)+3*sd(ukb_cauc_gene$genetic_heterozygosity.0.0)
rm_u_ids = ukb_cauc_gene$f.eid[which(ukb_cauc_gene$genetic_heterozygosity.0.0 >= uthr)]

# lower outliers
lthr = mean(ukb_cauc_gene$genetic_heterozygosity.0.0)-3*sd(ukb_cauc_gene$genetic_heterozygosity.0.0)
rm_l_ids = ukb_cauc_gene$f.eid[which(ukb_cauc_gene$genetic_heterozygosity.0.0 <= lthr)]

# both
rm_ids = sort(c(rm_u_ids, rm_l_ids))

ukb_cauc_gene = ukb_cauc_gene %>% filter(!f.eid %in% rm_ids)
dim(ukb_cauc_gene)

# write data
n_subs = as.character(nrow(ukb_cauc_gene))
out_file = paste0(base_dir, '/data/ukb/ukb_imaging_n',n_subs,'.csv')
write_csv(x=ukb_cauc_gene, path=out_file)

# write subjects for subsetting of genetic data
write.table(x=ukb_cauc_gene[c('f.eid','f.eid')], file=paste0(base_dir, '/data/ukb/genotyped/ukb_imaging_subs_before_cryptic_relatedness_n',n_subs,'.txt'), quote=F, row.names = F, col.names = F)


# ----------
# You should run genetic preprocessing in the "11_genetics_pipe.bash" script before running the code below
# ----------


# subset the imaging data to reflect non-cryptically related individuals
unrelated_subs = read.table(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n',n_subs,'_rm025.grm.id'))
ukb_cauc_gene_unrelated = ukb_cauc_gene[ukb_cauc_gene$f.eid %in% unrelated_subs$V1,]


# read eigenvalues and add to the data frame
eigenvecs = read.table(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.01_n',n_subs,'_rm025.eigenvec'))
colnames(eigenvecs) = c('FID','IID',paste0('eigen_maf01_',1:10))
ukb_cauc_gene_unrelated = merge(x=ukb_cauc_gene_unrelated, y=eigenvecs, by.x='f.eid', by.y='FID')


# read eigenvalues and add to the data frame
eigenvecs = read.table(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n',n_subs,'_rm025.eigenvec'))
colnames(eigenvecs) = c('FID','IID',paste0('eigen_maf05_',1:10))
ukb_cauc_gene_unrelated = merge(x=ukb_cauc_gene_unrelated, y=eigenvecs, by.x='f.eid', by.y='FID')


# write final csv file with imaging/phenotype/genetic data
new_n_subs = nrow(ukb_cauc_gene_unrelated)
write_csv(x=ukb_cauc_gene_unrelated, path=paste0(base_dir, '/data/ukb/ukb_cauc_gene_unrelated_n',new_n_subs,'.csv'))






