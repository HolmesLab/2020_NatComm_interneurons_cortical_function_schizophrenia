library(tidyverse)
library(psych)
library(doBy)
library(doParallel)

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path){
  schaff_dir       = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
  dscalar_template = paste0(schaff_dir, '/Schaefer2018_400Parcels_7Networks_order.dscalar.nii')
  parcel_info_file = paste0(schaff_dir, '/Schaefer2018_400Parcels_7Networks_order_info.txt')
  write_val_file   = paste0('/gpfs/milgram/project/holmes/kma52/tmp2.txt')
  write_delim(delim=' ', x=as.data.frame(values), path=write_val_file,col_names=F)
  save_path = out_path

  matfunc = 'plotVolOnSurface'
  cmd = paste0('/nexsan/apps/hpc/Apps/Matlab/R2014a/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/2018_interneuron/scripts\');',
                matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                write_val_file, '\',\'', save_path, '\'); exit;"')
  system(cmd)
}


# set up directories
proj_dir    = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ukb_mri_dir = '/gpfs/milgram/data/UKB/REPOSITORY/afni_fast'


# Atlas MNI ROI information
roi_labels = '/gpfs/milgram/project/holmes/HOLMES_UKB/atlas/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order.txt'
roi_table  = read.table(roi_labels)


# potential subjects
file_list = list.files(ukb_mri_dir)
subj_list = file_list[unlist(lapply(file_list, function(x) nchar(x)==7))]


# Read RSFA values from the Schaeffer 7-network atlas
rsfa_df           = as.data.frame( matrix(NA, nrow=length(subj_list), ncol=length(roi_table$V1)+1) )
colnames(rsfa_df) = c('UKB_ID', as.character(roi_table$V2)[1: length(roi_table$V2)])
i = 1
for ( subject in subj_list ){
  print(i)
  subj_dir  = paste0(ukb_mri_dir, '/', subject)
  rsfc_file = paste0(subj_dir, '/', subject, '_MNI_filtered_func_data_clean_despike_blur4_noGS_noWM_noCSF_Schaefer2018_400Parcels_7Net_FSLMNI152_2mm_RSFA.txt')

  # read RSFA values if the file exists
  if ( file.exists(rsfc_file) ){
    rsfa.dat = read_delim(rsfc_file, delim='\t')
    if ( length(t(rsfa.dat)) == 402 ){
      rsfa_df[i,] = c(subject, as.character(t(rsfa.dat)[3:402]))
    }
  } else {
    rsfa_df[i,1] = subject
  }
  i=i+1
}
output_file = paste0(proj_dir, '/data/ukb/rsfa/RSFA_MNI_filtered_func_data_clean_despike_blur4_noGS_noWM_noCSF_schaeff400.csv')
write.csv(x=rsfa_df, file=output_file, row.names = F)






# end














