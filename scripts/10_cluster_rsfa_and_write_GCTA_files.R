library(tidyverse)
library(factoextra)
library(caret)
library(mltools)
library(data.table)
library(Cairo)


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


# project directory
n = '9713'
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# read preprocessed combined pheno/rsfa data
csv_path = paste0(base_dir, '/data/ukb/ukb_cauc_gene_unrelated_n',as.character(n),'.csv')
ukb_subset = read_csv(csv_path)


# network parcel names
net_names = colnames(ukb_subset)[grep(paste0('RSFA_Net7'), colnames(ukb_subset))]
rsfa_only = ukb_subset[grep(paste0('RSFA_Net7'), colnames(ukb_subset))]


# regress each covariate from roi-based RSFA, save residuals
quant_covars = c('image_visit_age',
                 'BMI.2.0',
                 'height_standing',
                 'weight_preimaging.2.0', 
                 'lesion_vol',
                 'MRI_REST_motion.2.0',
                 'MRI_T1_invSNR.2.0',
                 'MRI_REST_invSNR_preproc.2.0',
                 'diastolic.2.0',
                 'systolic.2.0',
                 'MRI_T1_vetntCSF_norm.2.0', 
                 'MRI_T1_volGM_WM_headNorm.2.0',
                 'scanner_X_pos.2.0',
                 'scanner_Z_pos.2.0',
                 'scanner_Y_pos.2.0')

# scale quantitative covariates
scaled_qcovars = NULL
for (covar in quant_covars){
    new_covar = paste0(covar,'_scale')
    ukb_subset[[new_covar]] = as.numeric(scale(ukb_subset[[covar]]))
    scaled_qcovars = c(scaled_qcovars, new_covar)

    print(covar)
    print(summary(ukb_subset[[new_covar]]))
    print('')
}
ukb_subset$age2_scale = ukb_subset$image_visit_age_scale^2


# regress covariates from each parcel RSFA estimate
net_resid = NULL
for (net in paste0(net_names)){
  write(net,'')
  cur_resid = resid(lm(dat=ukb_subset, as.formula(paste0(net ,'~ age2_scale*sex_num.0.0 + image_visit_age_scale*sex_num.0.0 + sex_num.0.0 + UK_Biobank_assessment_centre.2.0 + age2_scale + ', paste(scaled_qcovars, collapse=' + ')))))
  net_resid = cbind(net_resid, cur_resid)
}
colnames(net_resid) = net_names
write.csv(x=net_resid, paste0(base_dir, '/data/ukb/rsfa/RSFA_n',as.character(n),'_residualized.csv'), quote=F, row.names=F)



# heirarchical clustering
set.seed(1234)
corr_dist  = dist(cor(net_resid))
clusters   = hclust(corr_dist, method = "ward.D2")
dend       = as.dendrogram(clusters)
clusterCut = cutree(clusters, 7)
cluster_df = data.frame(parcels=names(clusterCut), clusters=as.numeric(clusterCut))
write.csv(x=cluster_df, paste0(base_dir, '/data/ukb/rsfa/RSFA_n',as.character(n),'_cluster_7.csv'), quote=F, row.names=F)


# plot the clustering solution
i=7
clusterCut = cutree(clusters, i)

# set up paths/files for plotting
schaff_dir       = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/CBIG_private/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/HCP/fslr32k/cifti'
dscalar_template = paste0(schaff_dir, '/Schaefer2018_400Parcels_7Networks_order.dscalar.nii')
parcel_info_file = paste0(schaff_dir, '/Schaefer2018_400Parcels_7Networks_order_info.txt')
write_val_file   = paste0(base_dir, '/data/ukb/rsfa/clusters/RSFA_ward_clust', as.character(i),'_resid.txt')
write.table(x=as.numeric(clusterCut), file=write_val_file, row.names=F, col.names=F, quote=F)
save_path        = paste0(base_dir, '/data/ukb/rsfa/clusters/noscale_Schaff400_RSFA_ward_clust', as.character(i) ,'_n',as.character(n),'_resid.dlabel.nii')

# plot parcellation on surface
matfunc = 'plotVolOnSurface'
cmd = paste0('/nexsan/apps/hpc/Apps/Matlab/R2014a/bin/matlab -nosplash -nodesktop -r "cd(\'/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/scripts\');',
                matfunc, '(\'', dscalar_template, '\',\'', parcel_info_file, '\',\'',
                write_val_file, '\',\'', save_path, '\'); exit;"')
system(cmd) # plot

# add colors to the label file
label_path = paste0('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ukb/rsfa/clusters/cluster', as.character(i) ,'_colors.txt')
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_path, ' ', save_path))


# Plot between-subjects clustering of RSFA
pdf(width=5, height=1.5, file=paste0(base_dir, "/figures/PaperFig_RSFA_n",as.character(n),'_cluster7_dendrogram.pdf'))
clusters %>% fviz_dend(cex = 0.5, k = 7, palette = "jco")
dev.off()



# subjXroi matrix of func values
rsfa_only = ukb_subset[grep(paste0('RSFA_Net7'), colnames(ukb_subset))]


# average RSFA in each cortical cluster
clusterCut = cutree(clusters, 7)
avg_rsfa_noscale   = NULL
for (val in unique(clusterCut)){
  avg_rsfa_noscale = cbind(avg_rsfa_noscale, rowMeans(rsfa_only[clusterCut==val]))
}
avg_rsfa_noscale = as.data.frame(avg_rsfa_noscale)



# these names were determined from post-hoc visual inspection based on clustering profile in cortex
name_order = c('cingulo-opercular', 'visual', 'parietal', 'motor', 'prefrontal', 'limbicA','limbicB')
colnames(avg_rsfa_noscale) = name_order


# merge
ukb_subset_noscale = cbind(ukb_subset, avg_rsfa_noscale)


# write data files required for snp-wise heritability calculation (GCTA)
out_dir   = paste0(base_dir, '/data/ukb/rsfa')
pheno_dir = paste0(out_dir, '/pheno_RSFA')


# write cluster data for GCTA
cluster_RSFA = ukb_subset_noscale[c('f.eid','f.eid',name_order)]
cluster_RSFA_file = paste0(pheno_dir, '/RSFA_brit_n', as.character(n), '_cluster7_noscale.txt')
write.table(x=cluster_RSFA, file=cluster_RSFA_file, row.names=F, col.names=F, quote=F)

cluster_RSFA_colnames = paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_cluster7_colnames.txt')
write.table(x=name_order, file=cluster_RSFA_colnames, row.names=F, col.names=T, quote=F, )


# convert categorical covariates to binary
onehot_me     = ukb_subset[,c('sex_num.0.0','UK_Biobank_assessment_centre.2.0')]
onehot_me$sex_num.0.0 = as.factor(onehot_me$sex_num.0.0)
onehot_me$UK_Biobank_assessment_centre.2.0 = as.factor(onehot_me$UK_Biobank_assessment_centre.2.0)
one_hot_covar = as.data.frame(one_hot(data.table(onehot_me)))


# subjXroi matrix of func values
net_names = colnames(ukb_subset)[grep(paste0('RSFA_Net7'), colnames(ukb_subset))]
print(head(net_names))
pheno_dir = paste0(out_dir, '/pheno_RSFA')
if (dir.exists(pheno_dir) == F){
    dir.create(pheno_dir)
}


# write categorical covarites (GCTA)
covar_fields = c('f.eid','f.eid','sex_num.0.0','UK_Biobank_assessment_centre.2.0')
ccovar_path  = paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_ccovar.txt')
write.table(row.names=F, col.names=F, quote=F, x=ukb_subset[,covar_fields], file=ccovar_path)


# average RSFA in each individual cortical parcel (400)
write.table(x=ukb_subset[c('f.eid','f.eid',paste0(net_names))], file=paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_400parcel.txt'), row.names=F, col.names=F, quote=F)
write.table(x=net_names, file=paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_400parcel_colnames.txt'), row.names=F, col.names=F, quote=F)


# quantitative covariate file
regress_vars = c('f.eid','f.eid','image_visit_age_scale','BMI.2.0_scale','height_standing_scale','weight_preimaging.2.0_scale','MRI_T1_vetntCSF_norm.2.0_scale','MRI_T1_volGM_WM_headNorm.2.0_scale',
                  'MRI_REST_motion.2.0_scale','MRI_T1_invSNR.2.0_scale','MRI_REST_invSNR_preproc.2.0_scale', 'diastolic.2.0_scale', 'systolic.2.0_scale','age2_scale', 'lesion_vol_scale', 'scanner_X_pos.2.0_scale','scanner_Z_pos.2.0_scale','scanner_Y_pos.2.0_scale')
# z-transform each quantitative covar
qcovar = ukb_subset[regress_vars]
for (net in regress_vars[3:length(regress_vars)]){
  qcovar[net] = scale(qcovar[net])
}
# add eigenvectors
for (maf in c('01','05')){
    eigen_cols   = paste0('eigen_maf',maf,'_',1:10)
    write_qcovar = cbind(qcovar, ukb_subset[eigen_cols])
    write.table(x=write_qcovar, file=paste0(pheno_dir, '/RSFA_maf',maf,'_brit_n',as.character(n),'_qcovar.txt'), row.names=F, col.names=F, quote=F)

    # MAGMA covariate file
    write.table(x=cbind(qcovar, one_hot_covar, ukb_subset[eigen_cols]), file=paste0(pheno_dir, '/MAGMA_RSFA_maf_',maf,'_brit_n',as.character(n),'_covar.txt'), row.names=F, col.names=F, quote=F)
}
# average RSFA in each individual cortical parcel (400)
write.table(x=ukb_subset[c('f.eid','f.eid',paste0(net_names))], file=paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_400parcel.txt'), row.names=F, col.names=F, quote=F)
write.table(x=net_names, file=paste0(pheno_dir, '/RSFA_brit_n',as.character(n),'_400parcel_colnames.txt'), row.names=F, col.names=F, quote=F)



# Plot RSFA across cortex
nparcel='400'
net_num='7'
plot_matlab(values=as.numeric(colMeans(rsfa_only)), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_avg_RSFA_',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)


# parcel-wise expression of interneuron markers
schaeffer_mat = read_csv(paste0(base_dir, '/data/ahba/schaeffer400_7Net_expr_mat.csv'))
schaeff_sst   = schaeffer_mat[schaeffer_mat$gene == 'SST', 1:400]
schaeff_pvalb = schaeffer_mat[schaeffer_mat$gene == 'PVALB', 1:400]
schaeff_sst_pvalb = schaeff_sst-schaeff_pvalb


# read singlecell DFC data
nparcel = '400'
net_num = '7'
lake_dfc_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_',nparcel,'_',net_num,'Net_expr_mat.csv'))


# read singlecell VIS data
nparcel = '400'
net_num = '7'
lake_vis_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_VIS_',nparcel,'_',net_num,'Net_expr_mat.csv'))


# RSFA parcellation clusters
cluster_df = read_csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_n',as.character(n),'_cluster_7.csv'))


# combine parcel-wise expression and RSFA
parcel_array = clusterCut
expr_rsfa_df = NULL
RSFA_df      = as.data.frame(rsfa_only)
for (parcel in names(schaeff_sst)){
  write(parcel,'')
  expr_mean  = as.numeric(schaeff_sst_pvalb[parcel])
  pvalb_mean = as.numeric(schaeff_pvalb[parcel])
  sst_mean   = as.numeric(schaeff_sst[parcel])

  if ( is.na(expr_mean)){
    next
  } else {
      clust   = cluster_df$clusters[cluster_df$parcels %in% gsub('7Networks', 'RSFA_Net7', parcel)]
      cur_row = data.frame(sst_pvalb = expr_mean,
                          sst = sst_mean,
                          pvalb = pvalb_mean,
                          rsfa = mean(RSFA_df[[gsub('7Networks', 'RSFA_Net7', parcel)]]),
                          clust = as.numeric(clust),
                          name = parcel)
      dfc_cells    = as.data.frame(t(lake_dfc_dat[[parcel]]))
      colnames(dfc_cells) = lake_dfc_dat$gene
      cur_row      = cbind(cur_row, dfc_cells)
      vis_cells    = as.data.frame(t(lake_vis_dat[[parcel]]))
      colnames(vis_cells) = lake_vis_dat$gene
      cur_row      = cbind(cur_row, vis_cells)
      expr_rsfa_df = rbind(expr_rsfa_df, cur_row)
  }
}


# correlate RSFA with interneuron marker expression
cor.test(expr_rsfa_df$sst_pvalb, expr_rsfa_df$rsfa)
cor.test(expr_rsfa_df$sst_pvalb, expr_rsfa_df$rsfa, method='spearman')

cor.test(expr_rsfa_df$sst, expr_rsfa_df$rsfa)
cor.test(expr_rsfa_df$pvalb, expr_rsfa_df$rsfa)
cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$DFC_SST - expr_rsfa_df$DFC_PVALB)
cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$VIS_SST - expr_rsfa_df$VIS_PVALB)

cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$DFC_SST)
cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$VIS_SST)
cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$DFC_PVALB)
cor.test(expr_rsfa_df$rsfa, expr_rsfa_df$VIS_PVALB)



# all pairwise cell fraction subtractions correlated to RSFA
# Visual cortex
vis_cells   = colnames(expr_rsfa_df)[grep('VIS', colnames(expr_rsfa_df))]
vis_cells   = vis_cells[vis_cells != 'VIS_In3']
uniq_combos = combn(1:length(vis_cells),2)
vis_cell_diff_to_rsfa_df = NULL
for (combo in 1:ncol(uniq_combos)){
    # current cell combination
    cell_idx_a = uniq_combos[1,combo]
    cell_idx_b = uniq_combos[2,combo]
    
    # cell type names
    cell_a = vis_cells[cell_idx_a]
    cell_b = vis_cells[cell_idx_b]
    
    # subtract cell A expression from cell B
    cell_diff = scale(expr_rsfa_df[[cell_a]]) - scale(expr_rsfa_df[[cell_b]])
    cur_cor   = cor(expr_rsfa_df$rsfa, cell_diff)
    row_out   = data.frame(combo=paste0(cell_a, '_', cell_b), cor=cur_cor)
    vis_cell_diff_to_rsfa_df = rbind(vis_cell_diff_to_rsfa_df, row_out)
    
    # subtract cell B expression from cell A
    cell_diff = scale(expr_rsfa_df[[cell_b]]) - scale(expr_rsfa_df[[cell_a]])
    cur_cor   = cor(expr_rsfa_df$rsfa, cell_diff)
    row_out   = data.frame(combo=paste0(cell_b, '_', cell_a), cor=cur_cor)
    vis_cell_diff_to_rsfa_df = rbind(vis_cell_diff_to_rsfa_df, row_out)
}
# organize a bit
vis_plot_df     = as.data.frame(vis_cell_diff_to_rsfa_df)
vis_plot_df     = vis_plot_df[order(vis_plot_df$cor),]
vis_plot_df$idx = 1:nrow(vis_plot_df)

# plot distribution
vis_plot = ggplot(vis_plot_df, aes(idx, cor)) +
    geom_point(size=5, stroke=0, colour='black') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black",size=16),
          axis.text.y = element_text(colour="black",size=16),
          axis.title.x = element_text(colour="black",size=16),
          axis.title.y = element_text(colour="black",size=16)) +
    scale_y_continuous(expand = c(0, 0), limits=c(-.5,.5), breaks = seq(-.5, .5, by = .25))

library(CairoPDF)
vis_out = paste0(base_dir, '/figures/PaperFig_VIS_celldifferences_to_RSFA.pdf')
vis_out
CairoPDF(vis_out)
print(vis_plot)
dev.off()


# all pairwise cell fraction subtractions correlated to RSFA
# Frontal cortex
dfc_cells = colnames(expr_rsfa_df)[grep('DFC', colnames(expr_rsfa_df))]
uniq_combos = combn(1:length(dfc_cells),2)
dfc_cell_diff_to_rsfa_df = NULL
for (combo in 1:ncol(uniq_combos)){
    # current cell combination
    cell_idx_a = uniq_combos[1,combo]
    cell_idx_b = uniq_combos[2,combo]
    
    # cell type names
    cell_a = dfc_cells[cell_idx_a]
    cell_b = dfc_cells[cell_idx_b]
    
    # subtract cell A expression from cell B
    cell_diff = scale(expr_rsfa_df[[cell_a]]) - scale(expr_rsfa_df[[cell_b]])
    cur_cor   = cor(expr_rsfa_df$rsfa, cell_diff)
    row_out   = data.frame(combo=paste0(cell_a, '_', cell_b), cor=cur_cor)
    dfc_cell_diff_to_rsfa_df = rbind(dfc_cell_diff_to_rsfa_df, row_out)

    # subtract cell B expression from cell A
    cell_diff = scale(expr_rsfa_df[[cell_b]]) - scale(expr_rsfa_df[[cell_a]])
    cur_cor = cor(expr_rsfa_df$rsfa, cell_diff)
    row_out = data.frame(combo=paste0(cell_b, '_', cell_a), cor=cur_cor)
    dfc_cell_diff_to_rsfa_df = rbind(dfc_cell_diff_to_rsfa_df, row_out)
}
dfc_plot_df = as.data.frame(dfc_cell_diff_to_rsfa_df)
dfc_plot_df = dfc_plot_df[order(dfc_plot_df$cor),]
dfc_plot_df$idx = 1:nrow(dfc_plot_df)

# plot distribution
dfc_plot = ggplot(dfc_plot_df, aes(idx, cor)) +
    geom_point(size=5, stroke=0, colour='black') +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black",size=16),
          axis.text.y = element_text(colour="black",size=16),
          axis.title.x = element_text(colour="black",size=16),
          axis.title.y = element_text(colour="black",size=16)) +
    scale_y_continuous(expand = c(0, 0), limits=c(-.6,.6), breaks = seq(-.6, .6, by = .3))

dfc_out = paste0(base_dir, '/figures/PaperFig_DFC_celldifferences_to_RSFA.pdf')
dfc_out
CairoPDF(dfc_out)
print(dfc_plot)
dev.off()



# network colors
rgb_arr = c(rgb(42/255,204/255,160/255),
             rgb(124/255,18/255,133/255),
             rgb(230/255,147/255,36/255),
             rgb(70/255,130/255,182/255),
             rgb(205/255,61/255,84/255),
             rgb(219/255,248/255,167/255),
             rgb(122/255,135/255,50/255))


# plot expression RSFA to sst-pvalb
CairoPDF(width=4, height=4, file=paste0(base_dir, "/figures/PaperFig_schaeff400_clust7_RSFA_by_expression_plot.pdf"))
ggplot(expr_rsfa_df, aes(sst_pvalb, rsfa)) +
  geom_point(aes(colour = factor(clust)), size = 2, alpha=1) +
  scale_color_manual(values = rgb_arr) +
  theme_classic() +
  xlab('SST - PVALB expression') +
  ylab('RSFA') +
  expand_limits(x = c(-4, 6)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-4, 6, by = 2)) +
  expand_limits(y = c(0, 90)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 90, by = 15)) +
  geom_smooth(method = "lm", se = FALSE) + theme(legend.position = "none")
dev.off()



# Frontal cortex - each cell type fraction to RSFA
cells          = lake_dfc_dat$gene
dfc_cell_mat_t = as.data.frame(t(lake_dfc_dat[,1:400]))
colnames(dfc_cell_mat_t) = lake_dfc_dat$gene
dfc_cell_mat_t$region = gsub('7Networks', 'RSFA_Net7', rownames(dfc_cell_mat_t))

rsfa_means    = colMeans(rsfa_only)
rsfa_df       = data.frame(rsfa=rsfa_means, regions=names(rsfa_means))
cell_and_rsfa = merge(x=dfc_cell_mat_t, y=rsfa_df, by.x='region', by.y='regions')
cor_to_rsfa   = cor(cell_and_rsfa$rsfa, cell_and_rsfa[cells], use='pairwise.complete')

# organize for plotting
rsfa_cell_plot_df      = data.frame(cor2rsfa=as.numeric(cor_to_rsfa), cell=colnames(cor_to_rsfa))
rsfa_cell_plot_df      = rsfa_cell_plot_df[order(rsfa_cell_plot_df$cor2rsfa),]
rsfa_cell_plot_df$cell = factor(as.character(rsfa_cell_plot_df$cell), levels=as.character(rsfa_cell_plot_df$cell))
#
dfc_rsfa_cell_plot_df = rsfa_cell_plot_df
dfc_rsfa_cell_plot_df$cell = gsub('DFC_', '', dfc_rsfa_cell_plot_df$cell)
dfc_rsfa_cell_plot_df$ref  = 'DFC'

figure_out = paste0(base_dir,'/figures/PaperFig_celltype_DFC_cor_to_RSFA.pdf')
pdf(figure_out, width=1.5, height=2)
ggplot(data=rsfa_cell_plot_df, aes(x=cell, y=cor2rsfa)) +
    geom_bar(stat="identity") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_text(colour="black",size=6),
              axis.title.y = element_text(colour="black",size=6)) +
        scale_y_continuous(expand = c(0,.0), breaks=seq(-.5, .5, by=.25), limits=c(-.5,.5))
dev.off()




# Visual cortex - each cell type fraction to RSFA
cells          = lake_vis_dat$gene
vis_cell_mat_t = as.data.frame(t(lake_vis_dat[,1:400]))
colnames(vis_cell_mat_t) = lake_vis_dat$gene
vis_cell_mat_t$region = gsub('7Networks', 'RSFA_Net7', rownames(vis_cell_mat_t))

rsfa_means    = colMeans(rsfa_only)
rsfa_df       = data.frame(rsfa=rsfa_means, regions=names(rsfa_means))
cell_and_rsfa = merge(x=vis_cell_mat_t, y=rsfa_df, by.x='region', by.y='regions')
cor_to_rsfa   = cor(cell_and_rsfa$rsfa, cell_and_rsfa[cells], use='pairwise.complete')

# organize for plotting
rsfa_cell_plot_df      = data.frame(cor2rsfa=as.numeric(cor_to_rsfa), cell=colnames(cor_to_rsfa))
rsfa_cell_plot_df      = rsfa_cell_plot_df[order(rsfa_cell_plot_df$cor2rsfa),]
rsfa_cell_plot_df$cell = factor(as.character(rsfa_cell_plot_df$cell), levels=as.character(rsfa_cell_plot_df$cell))
rsfa_cell_plot_df      = rsfa_cell_plot_df[!is.na(rsfa_cell_plot_df$cor2rsfa),]

vis_rsfa_cell_plot_df = rsfa_cell_plot_df
vis_rsfa_cell_plot_df$cell = gsub('VIS_', '', vis_rsfa_cell_plot_df$cell)
vis_rsfa_cell_plot_df$ref  = 'VIS'

figure_out = paste0(base_dir,'/figures/PaperFig_celltype_VIS_cor_to_RSFA.pdf')
pdf(figure_out, width=1.5, height=2)
ggplot(data=rsfa_cell_plot_df, aes(x=cell, y=cor2rsfa)) +
    geom_bar(stat="identity") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_text(colour="black",size=6),
              axis.title.y = element_text(colour="black",size=6)) +
        scale_y_continuous(expand = c(0,.0), breaks=seq(-.5, .5, by=.25), limits=c(-.5,.5))
dev.off()



# combine the frontal/visual correlation values from above and plot both together
both_rsfa_cell_plot_df = rbind(dfc_rsfa_cell_plot_df, vis_rsfa_cell_plot_df)
both_rsfa_cell_plot_df = rbind(both_rsfa_cell_plot_df, data.frame(cor2rsfa=NA, cell='In2', ref='DFC'))
both_rsfa_cell_plot_df = rbind(both_rsfa_cell_plot_df, data.frame(cor2rsfa=NA, cell='In3', ref='VIS'))
both_rsfa_cell_plot_df = rbind(both_rsfa_cell_plot_df, data.frame(cor2rsfa=NA, cell='Ex2', ref='VIS'))
avg_cor = both_rsfa_cell_plot_df %>% group_by(cell) %>% summarise(x=mean(cor2rsfa, na.rm=T))
both_rsfa_cell_plot_df$cell = factor(both_rsfa_cell_plot_df$cell, levels=avg_cor$cell[order(avg_cor$x)])


figure_out = paste0(base_dir,'/figures/PaperFig_celltype_CELLTYPES_cor_to_RSFA.pdf')
figure_out
pdf(figure_out, width=2.5, height=1.5)
ggplot(data=both_rsfa_cell_plot_df, aes(x=cell, y=cor2rsfa, fill=ref)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
    #geom_bar(stat="identity", color="black", position=position_dodge(), bar_spacing=1) +
    #geom_errorbar(aes(ymin=norm_expr-expr_se, ymax=norm_expr+expr_se), width=.4, position=position_dodge(.9)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits=c(-.6,.6), breaks = seq(-.6, .6, by = .3)) +
    scale_fill_manual(values=c('#6E727F','#BEC1CC')) +
            theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=8, angle = 90, hjust = 1),
          axis.text.y = element_text(colour="black", size=8),
          axis.title.x = element_text(colour="black", size=8),
          axis.title.y = element_blank(),
            legend.position = "none")
dev.off()





# gene-based reference distributinos of pairwise gene subtractions to cortical RSFA

# read preprocessed UKB data
interneuron_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ukb_cauc_gene   = read_csv(paste0(interneuron_dir, '/data/ukb/ukb_cauc_gene_unrelated_n9713.csv'))
cluster_df      = read_csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_n',as.character(n),'_cluster_7.csv'))

# average parcelwise RSFA
net_names     = colnames(ukb_cauc_gene)[grep('RSFA', colnames(ukb_cauc_gene))]
avg_rsfa      = as.numeric(colMeans(ukb_cauc_gene[net_names]))
sst_min_pvalb = as.numeric(schaeff_sst) - as.numeric(schaeff_pvalb)

int_rsfa = cor.test(sst_min_pvalb, avg_rsfa)
uniq_combos = combn(1:nrow(schaeffer_mat),2)


# do lots of correlations - takes a bit
lower = 1
upper = 100000
gene_diff_to_rsfa = rep(NA, 1, ncol(uniq_combos))
while (lower < ncol(uniq_combos)){
    print(lower)
    if (upper > ncol(uniq_combos)){
        upper=ncol(uniq_combos)
    }
    expr_mat_a = schaeffer_mat[uniq_combos[1,lower:upper],1:400]
    expr_mat_b = schaeffer_mat[uniq_combos[2,lower:upper],1:400]
    chunk_cors = cor(t(expr_mat_a - expr_mat_b), avg_rsfa, use='pairwise.complete')
    gene_diff_to_rsfa[lower:upper] = chunk_cors
    lower=lower+100000
    upper=upper+100000
}
save(x=gene_diff_to_rsfa, file='/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/all_genepairs_w_rsfa.Rdata')
load('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/all_genepairs_w_rsfa.Rdata', verbose=T)


cor_df = data.frame(cors=gene_diff_to_rsfa)
length(which(cor_df$cors <= int_rsfa$estimate))
length(cor_df$cors)

pdf('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_SSTPVALB_diff_corr_to_RSFA_distribution.pdf',width=5,height=3)
ggplot(cor_df, aes(cors)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('All pairwise cortical spatial correlations')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=16),
          axis.title.y = element_text(colour="black", size=16)) +
    ylab("Density") +
    expand_limits(y=c(0,2)) +
    scale_x_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1,1,.5)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = int_rsfa$estimate, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()



# characterise gene pairs that are more associated to RSFA than SST-PVALB

# highly neg correlated gene combos
rsfa_neg_combos = uniq_combos[,which(cor_df$cors <= int_rsfa$estimate)]
x = table(schaeffer_mat$gene[c(rsfa_neg_combos[1,], rsfa_neg_combos[2,])])


# determine if the highly neg correlated gene pairs correlate with SST-PVALB
expr_mat_a = schaeffer_mat[rsfa_neg_combos[1,],1:400]
expr_mat_b = schaeffer_mat[rsfa_neg_combos[2,],1:400]
chunk_cors = cor(t(expr_mat_a - expr_mat_b), sst_min_pvalb, use='pairwise.complete')


# plot the distribution
df = data.frame(x=chunk_cors)
sst_plot = paste0(base_dir, '/figures/PaperFig_high_rsfa_genes_to_SST-PVALB.pdf')
sst_plot
pdf(sst_plot, height=2, width=3)
ggplot(df, aes(x)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('All pairwise cortical spatial correlations')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=16),
          axis.title.y = element_text(colour="black", size=16)) +
        scale_x_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1,1,.5))
dev.off()



# determine if genes in the highly RSFA related pairs correlate with SST and PVALB

# read gene ordering pvalb/sst
pvalb_df = read_csv(paste0(base_dir,'/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
sst_df   = read_csv(paste0(base_dir,'/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))

gradient_genes_to_pvalb = pvalb_df$pvalb_cor[pvalb_df$gene_symbol %in% names(rev(sort(x)))[1:500]]
gradient_genes_to_sst = sst_df$sst_cor[sst_df$gene_symbol %in% names(rev(sort(x)))[1:500]]
df = data.frame(x=gradient_genes_to_sst)

sst_plot = paste0(base_dir, '/figures/PaperFig_high_rsfa_genes_to_sst.pdf')
sst_plot
pdf(sst_plot, height=2, width=3)
ggplot(df, aes(x)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('All pairwise cortical spatial correlations')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=16),
          axis.title.y = element_text(colour="black", size=16)) +
        scale_x_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1,1,.5))
dev.off()

pvalb_plot = paste0(base_dir, '/figures/PaperFig_high_rsfa_genes_to_pvalb.pdf')
pvalb_plot
pdf(pvalb_plot, height=2, width=3)
df = data.frame(x=gradient_genes_to_pvalb)
ggplot(df, aes(x)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('All pairwise cortical spatial correlations')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=20),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=16),
          axis.text.y = element_text(colour="black", size=16),
          axis.title.x = element_text(colour="black", size=16),
          axis.title.y = element_text(colour="black", size=16)) +
        scale_x_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1,1,.5))
dev.off()



