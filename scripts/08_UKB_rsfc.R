library(tidyverse)
library(psych)
library(doBy)
library(doParallel)

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
    dscalar_template = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order.dscalar.nii')
    parcel_info_file = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcel_num,'Parcels_',net_num,'Networks_order_info.txt')
    write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/tmp/tmp.txt')
    write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
    save_path = out_path
    print(save_path)

    matfunc = 'plotVolOnSurface'
    cmd = paste0('/gpfs/milgram/apps/hpc.rhel7/software/MATLAB/2017b/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/scripts\');',
                matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                write_val_file, '\',\'', save_path, '\'); exit;"')
    system(cmd)
}



# set up directories
ukb_dir     = '/gpfs/milgram/data/UKB'
ukb_mri_dir = paste0(ukb_dir, '/REPOSITORY/afni_fast')
base_dir    = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
output_dir  = paste0(base_dir, '/data/ukb/rsfc')


# potential subjects
file_list = list.files(ukb_mri_dir)
subj_list = file_list[unlist(lapply(file_list, function(x) nchar(x)==7))]


# for each potential subject, try to read the connectivity matrix, and produce a flattened connectivity text file
mat_list = NULL
sub_arr  = NULL
ct = 1
cl = makeCluster(10)
registerDoParallel(cl)
foreach ( i = 1:length(subj_list) ) %dopar% {
  subject  = subj_list[i]
  subj_dir = paste0(ukb_mri_dir, '/', subject)
  out_file = paste0(base_dir, '/data/ukb/rsfc/roi430_anaticor/', subject, '_MNI_filtered_func_data_clean_despike_blur4_fanaticor_Parcel7net400CorticoThalStriatHippAmygSNVTA_noOverlap')
  if ( !file.exists(out_file)){
    rsfc_cort_file = paste0(subj_dir, '/', subject, '_MNI_filtered_func_data_clean_despike_blur4_fanaticor_Parcel7net400CorticoThalStriatHippAmygSNVTA_noOverlap_000.netcc')
    if ( file.exists(rsfc_cort_file) ){
      rsfc_table = read.table(rsfc_cort_file)
      rsfc_use   = rsfc_table[-c(1,2),]
      zscore     = rsfc_use[431:nrow(rsfc_use),]

      if ( nrow(zscore) == 430 & ncol(zscore) == 430 ){
        write.table(x=zscore[upper.tri(zscore)], file=out_file, row.names = F, quote=F, col.names = F)
        write(ct,'')
        sub_arr = c(sub_arr, subject)
        ct = ct + 1
      }
    }
  } else {
    write('file exists','')
  }
}
save(mat_list, file=paste0(output_dir, '/corticothalamicstriatal_roi_matrices_volrest.Rdata'))
save(sub_arr, file=paste0(output_dir, '/corticothalamicstriatal_roi_sub_arr_volrest.Rdata'))


# read the flat connectivity text files created in the previous processing loop
rsfc_dir = paste0(base_dir, '/data/ukb/rsfc/roi430_anaticor')
files    = list.files(rsfc_dir, pattern='_MNI_filtered_func_data_clean_despike_blur4_fanaticor_Parcel7net400CorticoThalStriatHippAmygSNVTA_noOverlap')
rsfc_mat = matrix(NA, nrow=length(files), ncol=92235)
sub_arr  = NULL
row = 1
for (f in files){
  write(row,'')
  subject = strsplit(f,'_')[[1]][1]
  sub_arr = c(sub_arr, subject)
  sub_dat = read.table(paste0(rsfc_dir, '/', f))
  rsfc_mat[row,] = as.numeric(sub_dat$V1)
  row = row + 1
}
rsfc_mat         = as_tibble(rsfc_mat)
rsfc_mat$sub_arr = sub_arr
save(rsfc_mat, file=paste0(base_dir, '/data/ukb/rsfc/roi430_anaticor_corticothalamicstriatal_roi_matrices_volrest.Rdata'))
load(paste0(base_dir, '/data/ukb/rsfc/roi430_anaticor_corticothalamicstriatal_roi_matrices_volrest.Rdata'), verbose=T)


# subset to include the same subjects used in the RSFA analyses
use_subs = read.table(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.fam'))
rsfc_mat = rsfc_mat[which(rsfc_mat$sub_arr %in% use_subs$V1),]


# plug the roiXroi correlation values back into a matrix
avg_edges = colMeans(rsfc_mat[paste0('V', 1:92235)])
avg_mat   = matrix(NA,430,430)
avg_mat[upper.tri(avg_mat)] = avg_edges
avg_mat   = avg_mat[1:430,1:430]

# mirror the corr matrix
avg_mat[lower.tri(avg_mat)] = t(avg_mat)[lower.tri(avg_mat)]


# add column names to the ROI marix
schaeffer_info = read.table('/gpfs/milgram/project/holmes/HOLMES_UKB/atlas/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order.txt')
striat_info    = c('striat_visual','striat_sommot','striat_dorsattn','striat_ventattn','striat_limbic','striat_control','striat_default')
thal_info      = c('thal_default','thal_cingulopercular','thal_motor','thal_control','thal_LatOcc',
                    'thal_medOcc','thal_medTemp','thal_temp','thal_dorsAttn')
crbl_info      = c('crbl_visual','crbl_sommot','crbl_dorsattn','crbl_ventattn','crbl_limbic','crbl_control','crbl_default')
col_names      = c(as.character(schaeffer_info$V2), striat_info, thal_info, crbl_info, c('lh_hippo','rh_hippo','lh_amyg','rh_amyg', 'SN_compacta','SN_retic','VTA'))


# add ROI names to average matrix
colnames(avg_mat) = col_names
rownames(avg_mat) = col_names
avg_mat = as.data.frame(avg_mat)


# write the parcelwise correlation map for each subcortical striatal/thalamus ROI
for ( net in c(striat_info, thal_info) ){
  cort_vals = avg_mat[[net]][1:400]
  plot_matlab(values=cort_vals, out_path=paste0(base_dir, '/figures/', net, '_anaticor_avg.dscalar.nii'))
}

# Create the SST dlabel.nii file
striat_limbic = avg_mat[['striat_limbic']][1:400]
striat_limbic[striat_limbic < 0.05]  = 0
striat_limbic[striat_limbic >= 0.05] = 2
thal_default <- avg_mat[['thal_default']][1:400]
thal_default[thal_default < 0.05]  = 0
thal_default[thal_default >= 0.05] = 1
limbic_composite = striat_limbic + thal_default # mPFC region with positive limbic striatal/thalamic ROI


ahba_parcel = read.csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/schaeffer400_7Net_expr_mat.csv')
sst_expr   = ahba_parcel[ahba_parcel$gene == 'SST',1:400]
pvalb_expr = ahba_parcel[ahba_parcel$gene == 'PVALB',1:400]
parcel_sst_pvalb = sst_expr - pvalb_expr

# plot the expression values for cortical ROIs that are positively correlated to both limbic striatum and default thalamus
schaeff_sst_pvalb = data.frame(sst_expr=as.numeric(t(sst_expr)), pvalb_expr=as.numeric(t(pvalb_expr)), parcel_sst_pvalb =as.numeric(t(parcel_sst_pvalb)), category = limbic_composite)
schaeff_sst_pvalb$category[schaeff_sst_pvalb$category != 3] = 1
plot_schaeff = schaeff_sst_pvalb[!is.na(schaeff_sst_pvalb$sst_expr),]

plot_summ    = summaryBy(parcel_sst_pvalb~category, plot_schaeff, FUN=summary)
plot_summ$category = as.factor(plot_summ$category)
summary(aov(parcel_sst_pvalb ~ category, data=plot_schaeff))


library(Cairo)
CairoPDF(paste0(base_dir, '/figures/sst_anaticor_striatothalamic_parcels_vs_all_others.pdf'), height=2, width=4)
ggplot(data=plot_summ, aes(x=category, y=parcel_sst_pvalb.Mean)) +
  geom_pointrange(aes(ymin=plot_summ[['parcel_sst_pvalb.1st Qu.']], ymax=plot_summ[['parcel_sst_pvalb.3rd Qu.']]), position=position_dodge(0.4), size=.8) +
  geom_linerange(aes(ymin=plot_summ[['parcel_sst_pvalb.Min.']], ymax=plot_summ[['parcel_sst_pvalb.Max.']]), position=position_dodge(0.4), size=.2) +
  theme_minimal() + scale_color_manual(values=c('#F15854','#5DA5DA')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20),
        axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
        axis.line.y = element_blank(),
        axis.text.y = element_text(colour="black",angle=90, hjust=1, size=5),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10)) +
  scale_y_continuous(expand = c(0,0), limits=c(-6,6)) +
  ylab("") +
  xlab("sst striato-thalamic expr vs all others") + coord_flip()
dev.off()


# striatal/thalamic limbic overlap dlabel file
plot_matlab(values=limbic_composite, out_path=paste0(base_dir, '/figures/striat_thal_limbic_anaticor.dlabel.nii'), parcel_num=400, net_num=7)
save_path  = paste0(base_dir, '/figures/striat_thal_limbic_anaticor.dlabel.nii')
label_path = paste0(base_dir, '/figures/limbic_striatothal_colors.txt')
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_path, ' ', save_path))



# Identify cortical regions with positive sommot striatal/thalamic connectivity
striat_sommot = avg_mat[['striat_sommot']][1:400]
striat_sommot[striat_sommot < 0.03]=0
striat_sommot[striat_sommot >= 0.03]=2
thal_motor = avg_mat[['thal_motor']][1:400]
thal_motor[thal_motor < 0.03]=0
thal_motor[thal_motor >= 0.03]=1
sommot_composite = striat_sommot + thal_motor

# plot the expression values for cortical ROIs that are positively correlated to both limbic striatum and default thalamus
schaeff_sst_pvalb = data.frame(sst_expr=as.numeric(t(sst_expr)), pvalb_expr=as.numeric(t(pvalb_expr)), parcel_sst_pvalb =as.numeric(t(parcel_sst_pvalb)), category = sommot_composite)
schaeff_sst_pvalb$category[schaeff_sst_pvalb$category != 3] = 1
plot_schaeff = schaeff_sst_pvalb[!is.na(schaeff_sst_pvalb$sst_expr),]

plot_summ    = summaryBy(parcel_sst_pvalb~category, plot_schaeff, FUN=summary)
plot_summ$category = as.factor(plot_summ$category)
summary(aov(parcel_sst_pvalb ~ category, data=plot_schaeff))



library(Cairo)
CairoPDF(paste0(base_dir, '/figures/pvalb_anaticor_striatothalamic_parcels_vs_all_others.pdf'), height=2, width=4)
ggplot(data=plot_summ, aes(x=category, y=parcel_sst_pvalb.Mean)) +
  geom_pointrange(aes(ymin=plot_summ[['parcel_sst_pvalb.1st Qu.']], ymax=plot_summ[['parcel_sst_pvalb.3rd Qu.']]), position=position_dodge(0.4), size=.8) +
  geom_linerange(aes(ymin=plot_summ[['parcel_sst_pvalb.Min.']], ymax=plot_summ[['parcel_sst_pvalb.Max.']]), position=position_dodge(0.4), size=.2) +
  theme_minimal() + scale_color_manual(values=c('#F15854','#5DA5DA')) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.y = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20),
        axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
        axis.line.y = element_blank(),
        axis.text.y = element_text(colour="black",angle=90, hjust=1, size=5),
        axis.text.x = element_text(colour="black",size=10),
        axis.title.y = element_text(colour="black",size=10),
        axis.title.x = element_text(colour="black",size=10)) +
  scale_y_continuous(expand = c(0,0), limits=c(-6,6)) +
  ylab("") +
  xlab("sst striato-thalamic expr vs all others") + coord_flip()
dev.off()



# striatal/thalamic sommot overlap dlabel file
plot_matlab(values=sommot_composite, out_path=paste0(base_dir, '/figures/striat_thal_sommot_anaticor.dlabel.nii'), parcel_num=400, net_num=7)
save_path  = paste0(base_dir, '/figures/striat_thal_limbic_anaticor.dlabel.nii')
label_path = paste0(base_dir, '/figures/limbic_striatothal_colors.txt')
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_path, ' ', save_path))




plot_matlab(values=sommot_composite, out_path=paste0(base_dir, '/figures/striat_thal_sommot_anaticor.dlabel.nii'))
save_path = paste0(base_dir, '/figures/striat_thal_sommot_anaticor.dlabel.nii')
label_path = paste0(base_dir, '/figures/sommot_striatothal_colors.txt')
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_path, ' ', save_path))


