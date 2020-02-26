library(tidyverse)
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


# Set up paths
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# bed/bim/fam file for computing PRS
target   = paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025')


# remove major histocompatability complex SNPs (except for top SCZ associated variant), because of complex LD in this area
scz_gwas = paste0(base_dir, '/data/scz_gwas/ckqny.scz2snpres')
mhc_file = paste0(base_dir, '/data/ukb/genotyped/MHC.range')
write.table(c("6 25000000 35000000 mhc"), mhc_file, col.names = F, row.names = F, quote = F)


# remove MHC, but keep top SNP
plink = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/plink_1.09/plink'


# read target SNP data
target_bim = fread(paste0(target,'.bim'))
target_bim = target_bim[order(target_bim$V2),]


# read SCZ GWAS data
base = fread(paste(scz_gwas,sep = "/"))


# identify MHC SNPs
MHC_component = base %>% filter(CHR==6 & BP >= 25000000 & BP <= 35000000)
MHC_component = MHC_component[which(MHC_component$SNP %in% target_bim$V2),]


# identify top-snp
top_snp    = MHC_component$SNP[which(MHC_component$P == min(MHC_component$P))]
remove_mhc = MHC_component$SNP[!MHC_component$SNP == top_snp]
write.table(remove_mhc, file = paste0(base_dir, "/data/ukb/genotyped/mhc_remove.txt"), col.names = F,row.names = F, quote = F)


# plink removal of MHC SNPs
plink_cmd = paste(plink," --noweb --bfile ", target," --exclude ",paste0(base_dir, "/data/ukb/genotyped/mhc_remove.txt")," --make-bed --out ", paste0(target, '_onlyTopMHC'), sep = "")
print(plink_cmd)
system(plink_cmd)


# Calculate PGRS with PRSice
sice_dir = '/gpfs/milgram/project/holmes/HOLMES_UKB/external/PRSice_v2.1.0'
bfile    = paste0(target, '_onlyTopMHC')


# Ripke Schizophrenia GWAS - formatted for PRSice
sice_scz_gwas = '/gpfs/milgram/project/holmes/HOLMES_UKB/projects/ukb_genomics/gwas_summ_stats/PRSice.PGC.scz2.Ripke.Nature2014.txt'


# create/run PRSice cmd
sice_cmd = paste(paste0(sice_dir, '/PRSice_linux'),
                   '--base', sice_scz_gwas,
                   '--target', bfile,
                   '--thread 4',
                   '--stat OR',
                   '--binary-target T',
                   '--all-score',
                   '--no-regress',
                   '--fastscore',
                   '--missing CENTER',
                   '--out', paste0(base_dir,'/data/prsice/scz_full_score_onlyTopMHC'))
print(sice_cmd)
system(sice_cmd)


# read polygenic risk score
gwas_df = read.table(paste0(base_dir,'/data/prsice/scz_full_score_onlyTopMHC.all.score'), header=T)
colnames(gwas_df) = c('FID','IID',paste0('prs_', c('0.001','0.05','0.1','0.2','0.3','0.4','0.5','1.0')))


# Schizophrenia PRS predicting cluster-wise RSFA
rsfa_in = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_cluster7_noscale.txt')
cluster = read_delim(rsfa_in, col_names = F, delim=' ')
colnames(cluster) = c('IID','FID', paste0('rsfa',1:7))


# Read UKB phenotype file
ukb_dat_in = paste0(base_dir, '/data/ukb/ukb_cauc_gene_unrelated_n9713.csv')
ukb_subset = read_csv(ukb_dat_in)
colnames(ukb_subset) = gsub('7Net', 'net7Net', colnames(ukb_subset)) # columns shouldnt start with numbers


quant_covars = c('image_visit_age','BMI.2.0','height_standing','weight_preimaging.2.0', 'lesion_vol',
                  'MRI_REST_motion.2.0','MRI_T1_invSNR.2.0','MRI_REST_invSNR_preproc.2.0', 'diastolic.2.0', 'systolic.2.0',
                  'MRI_T1_vetntCSF_norm.2.0', 'MRI_T1_volGM_WM_headNorm.2.0',
                  'scanner_X_pos.2.0','scanner_Z_pos.2.0','scanner_Y_pos.2.0')
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
ukb_subset$age2 = ukb_subset$image_visit_age^2


# merge PRS with UKB phenotype
rsfa_gwas = merge(ukb_subset, gwas_df, by.x='f.eid', by.y='FID')


# merge PRS/UKB phenotype with RSFA clusters
rsfa_gwas = merge(rsfa_gwas, cluster, by.x='f.eid', by.y='IID')


# check that we're examining the correct subjects
use_subs  = read.table(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.grm.id'), header=F)
rsfa_gwas = rsfa_gwas[which(rsfa_gwas$FID %in% use_subs$V1),]
print(dim(rsfa_gwas))


# SCZ-PRS predicting cluster-wise RSFA
prs_mat    = NULL
net_arr    = c('cingulo-opercular', 'visual', 'parietal', 'motor', 'prefrontal', 'limbicA', 'limbicB')
thresh_arr = c('prs_1.0','prs_0.5','prs_0.1','prs_0.05')
PCs        = paste0('eigen_maf05_', 1:10, collapse=' + ')
for (thresh in thresh_arr){
  for ( i in 1:7 ){

    covariates = paste0('age2_scale*sex_num.0.0 + image_visit_age_scale*sex_num.0.0 + sex_num.0.0 + UK_Biobank_assessment_centre.2.0 + age2_scale + ', paste(scaled_qcovars, collapse=' + '), '+', PCs)
    out        = summary(lm(formula( paste0('scale(rsfa',i,') ~ scale(',thresh,') + ', covariates)), rsfa_gwas))
    prs_vals   = as_tibble(t(out$coefficients[rownames(out$coefficients) == paste0('scale(',thresh,')'),]))
    colnames(prs_vals)  = c('est','se','tval','pval')
    prs_vals$pgs_thresh = thresh
    prs_vals$net = net_arr[i]
    prs_mat      = rbind(prs_mat, prs_vals)
  }
}
prs_mat %>% filter(pgs_thresh == 'prs_1.0')
write.csv(x=prs_mat, paste0(base_dir, '/figures/PaperData_stats_RSFAcluster_scz_prs.csv'), row.names=F, quote=F)


# print FDR-BH corrected stats for each PRS-pval threshold
for (thresh in thresh_arr){
    print(thresh)
    top_cors      = prs_mat[prs_mat$pgs_thresh == thresh,]
    top_cors$qval = p.adjust(top_cors$pval, method='BH')
    print(top_cors[top_cors$qval < .05,])
    print(top_cors[top_cors$qval < .1,])
}


# plot PRS-RSFA relationship
prs_mat$net        = factor(prs_mat$net, levels=c('limbicA', 'limbicB', 'cingulo-opercular', 'parietal', 'prefrontal', 'motor', 'visual'))
prs_mat$pgs_thresh = factor(prs_mat$pgs_thresh, levels=c('prs_1.0','prs_0.5','prs_0.1','prs_0.05'))
fig_out = paste0(base_dir, '/figures/PaperFig_parcel7_scz_gwas_rsfa_new.pdf')
fig_out

pdf(fig_out, width=4, height=2.25)
ggplot(prs_mat, aes(x=net, y=est, fill = pgs_thresh)) +
  geom_bar(aes(fill = pgs_thresh), stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(fill = pgs_thresh, ymin=est-se, ymax=est+se), width=.5,  position=position_dodge(.9)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(-.04,.02, .01), expand = c(0, 0), limits=c(-.04,.02))
dev.off()



# Parcel-wise PRS score to RSFA analysis
parcel_prs_mat = NULL
net_names = colnames(rsfa_gwas)[grep('RSFA_', colnames(rsfa_gwas))]
thresh    = 'prs_1.0'
PCs       = paste0('eigen_maf05_', 1:10, collapse=' + ')
for ( net in net_names ){
    write(net, '')
    covariates = paste0('age2_scale*sex_num.0.0 + image_visit_age_scale*sex_num.0.0 + sex_num.0.0 + UK_Biobank_assessment_centre.2.0 + age2_scale + ', paste(scaled_qcovars, collapse=' + '), '+', PCs)
    out        = summary(lm(formula( paste0('scale(',net,') ~ scale(',thresh,') + ', covariates)), rsfa_gwas))

    prs_vals = as_tibble(t(out$coefficients[rownames(out$coefficients) == paste0('scale(',thresh,')'),]))
    colnames(prs_vals)  = c('est','se','tval','pval')
    prs_vals$pgs_thresh = thresh
    prs_vals$net   = net
    parcel_prs_mat = rbind(parcel_prs_mat, prs_vals)
}

# combine regression DF with info about parcel-wise RSFA heritability
parcel_prs_mat$idx = 1:nrow(parcel_prs_mat)
hsq_overall    = read_csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv')
parcel_prs_mat = merge(x=parcel_prs_mat, y=hsq_overall, by='net')
parcel_prs_mat = parcel_prs_mat[order(parcel_prs_mat$idx.x),]

parcel_prs_mat$h2snp = parcel_prs_mat$tot_hsq

# read/combind dat on Parcelwise SCZ effects to gene expression
schaeff_mat = read_csv(paste0(base_dir, '/data/ahba/schaeffer400_7Net_expr_mat.csv'))


pvalb_exr = t(schaeff_mat[schaeff_mat$gene == 'PVALB',1:400])
rownames(pvalb_exr) = gsub('7Networks_', 'RSFA_Net7_', rownames(pvalb_exr))

colnames(pvalb_exr) == parcel_prs_mat$net
cor.test(as.numeric(pvalb_exr), parcel_prs_mat$est)

# Plot parcel-wise SCZ PRS-RSFA relationship
save_path = paste0(base_dir, '/figures/PaperFig_SCZ_GWAS_prs1_0_schaeff400_parcel_wise_corr.dscalar.nii')
plot_matlab(values=as.numeric(parcel_prs_mat$est), out_path=save_path, parcel_num='400', net_num='7')



# get rid of NA expression values
schaeff_nonan = schaeff_mat[,which(!is.na(schaeff_mat[1,1:400]))]
est_nonan     = parcel_prs_mat$est[which(!is.na(schaeff_mat[1,1:400]))]

# gene-wise correlation to SCZ-RSFA effects
cor_arr      = cor(est_nonan, t(schaeff_nonan))
prs_ahba_cor = data.frame(cor=as.numeric(cor_arr), gene=schaeff_mat$gene)
prs_ahba_cor_ordered = prs_ahba_cor[order(prs_ahba_cor$cor),]
which(prs_ahba_cor_ordered$gene == 'PVALB')/length(prs_ahba_cor_ordered$gene)
which(prs_ahba_cor_ordered$gene == 'PVALB')

# same as above, but with spearman
sp_cor_arr      = cor(est_nonan, t(schaeff_nonan), method='spearman')
sp_prs_ahba_cor = data.frame(cor=as.numeric(sp_cor_arr), gene=schaeff_mat$gene)
sp_prs_ahba_cor_ordered = sp_prs_ahba_cor[order(sp_prs_ahba_cor$cor),]
which(sp_prs_ahba_cor_ordered$gene == 'PVALB')/length(sp_prs_ahba_cor_ordered$gene)
which(sp_prs_ahba_cor_ordered$gene == 'PVALB')


# Extract pvalb/sst expression
# pvalb
pvalb_expr = as.numeric(schaeff_mat[schaeff_mat$gene == 'PVALB',1:400])
pvalb_expr[is.na(pvalb_expr)] = 0
# sst
sst_expr =  as.numeric(schaeff_mat[schaeff_mat$gene == 'SST',1:400])
sst_expr[is.na(sst_expr)] = 0


# cluster assignment for each parcel (for coloring points in plot)
parcels = read.csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_n9713_cluster_7.csv'))

# put all the expression/prs/cluster data together for plotting
plot_parcel_prs_mat            = parcel_prs_mat
plot_parcel_prs_mat$pvalb_expr = pvalb_expr
plot_parcel_prs_mat$sst_expr   = sst_expr
plot_parcel_prs_mat$cluster    = parcels$clusters

# get rid of missing data
plot_parcel_prs_mat = plot_parcel_prs_mat[plot_parcel_prs_mat$pvalb_expr != 0,]

# pvalb predicting scz-effect, controlling for heritability
cor.test(plot_parcel_prs_mat$pvalb_expr, plot_parcel_prs_mat$est)
cor.test(plot_parcel_prs_mat$pvalb_expr, plot_parcel_prs_mat$est, method='spearman')
out = summary(lm('scale(est) ~ scale(pvalb_expr) + scale(h2snp)', plot_parcel_prs_mat))
out

# color information for hsq/expression plots
net_arr = c('cingulo-opercular', 'visual', 'parietal', 'motor', 'prefrontal', 'limbicA', 'limbicB')
rgb_arr = c(rgb(42/255,204/255,160/255),
             rgb(124/255,18/255,133/255),
             rgb(230/255,147/255,36/255),
             rgb(70/255,130/255,182/255),
             rgb(205/255,61/255,84/255),
             rgb(122/255,135/255,50/255),
             rgb(219/255,248/255,167/255))


# create SCZ-PRS/expression plots
for ( plot_vals in c('pvalb_expr','sst_expr') ){

  # get rid of extreme heritability values (noise)
  parcel_hsq_subset_noextrema = plot_parcel_prs_mat

  # extract the current expression data (i.e. pvalb/sst) to examine
  parcel_hsq_subset_noextrema$plot_expr = parcel_hsq_subset_noextrema[[plot_vals]]

  # remove expression outliers
  upper = mean(parcel_hsq_subset_noextrema$plot_expr) + 4*sd(parcel_hsq_subset_noextrema$plot_expr)
  lower =  mean(parcel_hsq_subset_noextrema$plot_expr) - 4*sd(parcel_hsq_subset_noextrema$plot_expr)
  parcel_hsq_subset_noextrema = parcel_hsq_subset_noextrema %>% filter(plot_expr <= upper & plot_expr >= lower)

  # print correlation between herit/expression
  print(cor.test(parcel_hsq_subset_noextrema$plot_expr, parcel_hsq_subset_noextrema$est ))

  xrange = c(-3,3)
  if (plot_vals == 'pvalb_expr'){
    xrange=c(-3, 3)
  }
    
  # plot expression by heritability
  out_file = paste0(base_dir, '/figures/PaperFig_scz_prs_by_',plot_vals,'_plot.pdf')
  write(out_file,'')
  CairoPDF(width=2.5, height=2.5, file=out_file)
  p = ggplot(parcel_hsq_subset_noextrema[c('plot_expr','est','cluster')], aes(plot_expr, est)) +
    geom_point(aes(colour = factor(cluster)), size = 1.5, alpha=1) +
    scale_color_manual(values = rgb_arr) +
    theme_classic() +
    xlab(plot_vals) +
    ylab(paste0('prs^2')) +
    expand_limits(x = xrange) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(xrange[1], xrange[2], by = 1)) +
    expand_limits(y = c(-0.04, .02)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-0.04, .02, by = .02)) +
    geom_smooth(method = "lm", se = FALSE) +
    theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), legend.position = "none")
  print(p)
  dev.off()
}


pvalb_cor = prs_ahba_cor_ordered$cor[prs_ahba_cor_ordered$gene == 'PVALB']
mean_cor = mean(prs_ahba_cor_ordered$cor)

CairoPDF('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/scz_prs_map_to_all_genes.pdf',family='Arial',width=5,height=3)
ggplot(prs_ahba_cor_ordered, aes(cor)) +
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
    geom_vline(xintercept = pvalb_cor, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = mean_cor, color='black', size=1, linetype='longdash')
dev.off()



# Correspondence bewteen Single cell deconvolution maps and SCZ-RSFA effects

# DFC - create averaged plots
nparcel = '400'
net_num = '7'
lake_dfc_dat   = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_',nparcel,'_',net_num,'Net_expr_mat.csv'))
lake_dfc_dat_t = as.data.frame(t(lake_dfc_dat[,1:400]))

# row/col names
colnames(lake_dfc_dat_t) = lake_dfc_dat$gene
lake_dfc_dat_t$net = rownames(lake_dfc_dat_t)
lake_dfc_dat_t$net = gsub('7Networks','RSFA_Net7', lake_dfc_dat_t$net)

# correlate cell deconvolution to SCZ PRS
DFC_parcel_prs_mat = merge(x=parcel_prs_mat, y=lake_dfc_dat_t, by='net')
dfc_singlecell_cors_to_PRS = cor(DFC_parcel_prs_mat$est, DFC_parcel_prs_mat[lake_dfc_dat$gene], use='pairwise.complete')
dfc_cell_df = data.frame(cor2prs=as.numeric(t(dfc_singlecell_cors_to_PRS)), cell=rownames(t(dfc_singlecell_cors_to_PRS)))
dfc_cell_df$cell = gsub('DFC_', '', dfc_cell_df$cell)

# spatial correlation of SCZ-PRS effects to DFC derived cell esimates
cor.test(DFC_parcel_prs_mat$est, DFC_parcel_prs_mat$DFC_PVALB, use='pairwise.complete')
cor.test(DFC_parcel_prs_mat$est, DFC_parcel_prs_mat$DFC_PVALB, use='pairwise.complete', method='spearman')



# VISUAL cell reference
lake_vis_dat   = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_VIS_',nparcel,'_',net_num,'Net_expr_mat.csv'))
lake_vis_dat_t = as.data.frame(t(lake_vis_dat[,1:400]))

# row/col names
colnames(lake_vis_dat_t) = lake_vis_dat$gene
lake_vis_dat_t$net = rownames(lake_vis_dat_t)
lake_vis_dat_t$net = gsub('7Networks','RSFA_Net7', lake_vis_dat_t$net)
lake_vis_dat_t$VIS_In3 = NULL

# correlate cell deconvolution to SCZ PRS
VIS_parcel_prs_mat = merge(x=parcel_prs_mat, y=lake_vis_dat_t, by='net')
cells = colnames(lake_vis_dat_t)[1:(ncol(lake_vis_dat_t)-1)]
vis_singlecell_cors_to_PRS = cor(VIS_parcel_prs_mat$est, VIS_parcel_prs_mat[cells], use='pairwise.complete')
vis_singlecell_cors_to_PRS = as.data.frame(vis_singlecell_cors_to_PRS)

# spatial correlation of SCZ-PRS effects to VIS derived cell esimates
cor.test(VIS_parcel_prs_mat$est, VIS_parcel_prs_mat$VIS_PVALB, use='pairwise.complete')
cor.test(VIS_parcel_prs_mat$est, VIS_parcel_prs_mat$VIS_PVALB, use='pairwise.complete', method='spearman')


# spatial correlation of each cells spatial map to SCZ-PRS effects
vis_cell_df = data.frame(cor2prs=as.numeric(t(vis_singlecell_cors_to_PRS)), cell=rownames(t(vis_singlecell_cors_to_PRS)))
vis_cell_df$cell = gsub('VIS_', '', vis_cell_df$cell)

dfc_cell_df$ref = 'DFC'
vis_cell_df$ref = 'VIS'
both_cell_df = rbind(dfc_cell_df, vis_cell_df)
both_cell_df = rbind(both_cell_df, data.frame(cor2prs=NA,cell='In2',ref='DFC'))
both_cell_df = rbind(both_cell_df, data.frame(cor2prs=NA,cell='Ex2',ref='VIS'))
both_cell_df = rbind(both_cell_df, data.frame(cor2prs=NA,cell='In3',ref='VIS'))


# plot in order of average spatial correlatoni to SCZ-PRS RSFA effects
avg_prs_cors = both_cell_df %>% group_by(cell) %>% summarise(mean_cor=mean(cor2prs, na.rm=T))
cell_order   = avg_prs_cors$cell[order(avg_prs_cors$mean_cor)]
both_cell_df$cell = factor(both_cell_df$cell, levels=cell_order)


figure_out = paste0(base_dir,'/figures/PaperFig_celltypes_to_scz_prs_map.pdf')
figure_out
pdf(figure_out, width=3, height=1.5)
ggplot(data=both_cell_df, aes(x=cell, y=cor2prs, fill=ref)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
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




