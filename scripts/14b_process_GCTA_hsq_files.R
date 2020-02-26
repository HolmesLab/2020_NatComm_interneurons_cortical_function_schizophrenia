library(tidyverse)
library(Cairo)
library(RColorBrewer)

# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir         = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
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
base_dir     = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
part_dir     = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/hsq_partition')
analysis_dir = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA')
parcel_names = read_csv(paste0(analysis_dir, '/RSFA_brit_n9713_400parcel_colnames.txt'), col_names=FALSE)


# Read heritability estimates of each Cluster
parcel_hsq   = NULL
analysis_dir = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/hsq')
for (region in c('cingulo_opercular','visual','parietal','motor','limbicB','prefrontal', 'limbicA')){
    hsq_file = paste0(analysis_dir, '/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_',region,'_clust7_noscale.hsq')
    hsq_dat  = read_delim(hsq_file, delim='\t')

    hsq  = hsq_dat[hsq_dat$Source == 'V(G)/Vp',]$Variance
    se   = hsq_dat[hsq_dat$Source == 'V(G)/Vp',]$SE
    pval = hsq_dat[hsq_dat$Source == 'Pval',]$Variance

    row = data.frame(region=region, maf='0.05', hsq=hsq, se=se, pval=pval, func_type = 'RSFA')
    parcel_hsq = rbind(parcel_hsq, row)
}



# order the networks for plotting
positions = c('limbicA','limbicB','cingulo_opercular','parietal','prefrontal','motor','visual')
color_arr = c('#DBF8A7','#7A8732','#2ACCA0','#E69324','#CD3D54','#4682B6','#7B297F')
parcel_hsq$color = factor(parcel_hsq$region, levels=positions)
plot_me = parcel_hsq

pdf(width=4, height=2, file=paste0(base_dir, "/figures/PaperFig_overall_clust7_allfunctypes_h2snp_partition.pdf"))
ggplot(plot_me, aes(x=region, y=hsq, fill=color)) +
  geom_bar( stat="identity",
            position=position_dodge()) +
  geom_errorbar(aes(ymin=hsq-se, ymax=hsq+se), width=.4,
                position=position_dodge(.9)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(breaks = seq(0,.4,.1), expand = c(0, 0), limits=c(0,.4)) +
  scale_fill_manual(values = color_arr) +
  scale_x_discrete(limits = positions)
dev.off()



# Parcel-wise heritability of RSFA
analysis_dir = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/hsq')
schaeffer_parcels = read.table('/gpfs/milgram/project/holmes/HOLMES_UKB/atlas/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order.txt')
roi_hsq = NULL
for (iter_parcel in schaeffer_parcels$V2){
    parcel   = gsub('7Networks','Net7',iter_parcel)
    hsq_file = paste0(analysis_dir, '/RSFA_',parcel,'_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.hsq')
    hsq_dat  = read_delim(hsq_file, delim='\t')

    hsq = hsq_dat[hsq_dat$Source == 'V(G)/Vp',]$Variance
    se  = hsq_dat[hsq_dat$Source == 'V(G)/Vp',]$SE

    row = data.frame(region=parcel, maf='05', hsq=hsq, se=se, func_type = 'RSFA')
    roi_hsq = rbind(roi_hsq, row)
}

# plot parcelwise RSFA on surface
nparcel = '400'
net_num = '7'
plot_matlab(values=as.numeric(roi_hsq$hsq), out_path=paste0(base_dir, '/figures/surface_plots/PaperPlot_rsfa_overall_hsq.dscalar.nii'), parcel_num=nparcel, net_num=net_num)



# read partitioned heritability for each of seven clusters
parcel_part_hsq = NULL
analysis_dir = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/hsq_partition')
for (region in c('cingulo-opercular','visual','parietal','motor','limbicB','prefrontal', 'limbicA')){
    hsq_file = paste0(analysis_dir, '/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_',region,'_clust7_noscale.hsq')

    hsq_dat   = read_delim(hsq_file, delim='\t')
    pvalb_hsq = hsq_dat[hsq_dat$Source == 'V(G1)/Vp',]$Variance
    sst_hsq   = hsq_dat[hsq_dat$Source == 'V(G2)/Vp',]$Variance
    pvalb_se  = hsq_dat[hsq_dat$Source == 'V(G1)/Vp',]$SE
    sst_se    = hsq_dat[hsq_dat$Source == 'V(G2)/Vp',]$SE
    tot_hsq   = hsq_dat[hsq_dat$Source == 'Sum of V(G)/Vp',]$Variance
    row = data.frame(region=region,
                      pvalb_hsq=pvalb_hsq,
                      sst_hsq=sst_hsq,sst_se=sst_se,
                      tot_hsq=tot_hsq,pvalb_se=pvalb_se,
                      norm_pvalb_hsq=pvalb_hsq/tot_hsq,
                      norm_sst_hsq=sst_hsq/tot_hsq, func_type='RSFA')
    parcel_part_hsq = rbind(parcel_part_hsq, row)
}
# organize partial hsq dataframe
pvalb = parcel_part_hsq[c('region','norm_pvalb_hsq','pvalb_se', 'pvalb_hsq')]
colnames(pvalb) = c('region', 'norm_hsq','partial_se','raw_hsq')
pvalb$gene = 'pvalb'

sst   = parcel_part_hsq[c('region','norm_sst_hsq','sst_se','sst_hsq')]
colnames(sst) = c('region', 'norm_hsq','partial_se','raw_hsq')
sst$gene = 'sst'

# combine into DF
part_hsq_plot = rbind(pvalb, sst)
head(part_hsq_plot)


# Calculate enrichment and p-values for each partition
part_hsq_plot$g_frac = 0
pvalb_snps = 9571 # these SNP counts were taken directly from the GCTA GRM calculation
sst_snps   = 8308 # these SNP counts were taken directly from the GCTA GRM calculation
total_snps = 337501

# fraction of sst/pvalb SNPs, relative to rest of SNP set
part_hsq_plot$g_frac[part_hsq_plot$gene=='pvalb'] = pvalb_snps/total_snps
part_hsq_plot$g_frac[part_hsq_plot$gene=='sst']   = sst_snps/total_snps

# normalized heritability, divided by genome fraction
part_hsq_plot$enrichment_fraction = part_hsq_plot$norm_hsq/part_hsq_plot$g_frac
part_hsq_plot$enrichment_se       = part_hsq_plot$partial_se/part_hsq_plot$g_frac


# calculate p-values for each partition
part_hsq_plot$zstat = part_hsq_plot$raw_hsq / part_hsq_plot$partial_se
part_hsq_plot$enrich_pval = 1-(0.5+0.5*pchisq(part_hsq_plot$zstat^2,1))
pvalb_enrich      = part_hsq_plot %>% filter(gene=='pvalb')
pvalb_enrich$qval = p.adjust(pvalb_enrich$enrich_pval, method='BH')

# print results
pvalb_enrich[pvalb_enrich$qval < .05,]
pvalb_enrich$enrichment_fraction - pvalb_enrich$enrichment_se*1.96

# SST
sst_enrich      = part_hsq_plot %>% filter(gene=='sst')
sst_enrich$qval = p.adjust(sst_enrich$enrich_pval, method='BH')

# print results
sst_enrich[sst_enrich$qval < .05,]
sst_enrich$enrichment_fraction - sst_enrich$enrichment_se*1.96



# order the networks for plotting
positions = c('limbicA','limbicB','cingulo-opercular','prefrontal','parietal','motor','visual')
part_hsq_plot$region = factor(part_hsq_plot$region, levels=positions)
color_arr = c('#DBF8A7','#7A8732','#2ACCA0','#CD3D54','#E69324','#4682B6','#7C8583')


# FIgure 5b
fig_out = paste0(base_dir, '/figures/PaperFig_clust7_noscale_partitioned_hsq_interneuron.pdf')
fig_out
pdf(fig_out, width=5, height=3)
p1 = ggplot(part_hsq_plot, aes(x=region, y=norm_hsq)) +
  geom_bar(aes(fill = gene), stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(fill = gene, ymin=norm_hsq-partial_se, ymax=norm_hsq+partial_se), width=.4,  position=position_dodge(.9)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = c( 0, .05, .1, .15, .20, .25), expand = c(0, 0), limits=c(-0.01,.25)) +
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
print(p1)
dev.off()




# Read partitioned heritability for each individual parcel
genes        = c('pvalb','sst')
roi_part_hsq = NULL
part_dir     = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/hsq_partition')
parcel_names = read_csv(paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_400parcel_colnames.txt'), col_names=FALSE)

for (i in 1:400){
  write(i,'')
  parcel = parcel_names$X1[i]
  hsq_in = paste0(part_dir, '/v2_',parcel, '_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025_1_500.hsq')

  if (!file.exists(hsq_in)){
    row = data.frame(net=as.character(parcel),
                      idx=i,
                      pvalb_hsq=0,
                      sst_hsq=0,
                      tot_hsq=0,
                      norm_pvalb_hsq=0,
                      norm_sst_hsq=0, func_type='RSFA')
  } else {
    hsq_dat   = read_delim(hsq_in, delim='\t')
    pvalb_hsq = hsq_dat[hsq_dat$Source == 'V(G1)/Vp',]$Variance
    sst_hsq   = hsq_dat[hsq_dat$Source == 'V(G2)/Vp',]$Variance
    tot_hsq   = hsq_dat[hsq_dat$Source == 'Sum of V(G)/Vp',]$Variance
    row = data.frame(net=as.character(parcel),
                      idx=i,
                      pvalb_hsq=pvalb_hsq,
                      sst_hsq=sst_hsq,
                      tot_hsq=tot_hsq,
                      norm_pvalb_hsq=pvalb_hsq/tot_hsq,
                      norm_sst_hsq=sst_hsq/tot_hsq,func_type='RSFA')
  }
    roi_part_hsq = rbind(roi_part_hsq, row)
}
write_csv(x=roi_part_hsq, path=paste0(base_dir, '/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv'))
roi_part_hsq = read_csv(paste0(base_dir, '/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv'))

# read schaeffer averaged expression data
schaeff_mat  = read_csv(paste0(base_dir, '/data/ahba/schaeffer400_7Net_expr_mat.csv'))
expr_parcels = colnames(schaeff_mat)[1:400]
gsub('7Networks_', 'RSFA_Net7_', expr_parcels) == roi_part_hsq$net


# PVALB expression across each Scheafer parcel
pvalb_expr = as.numeric(schaeff_mat[schaeff_mat$gene == 'PVALB',1:400])
pvalb_expr[is.na(pvalb_expr)] = 0

# SST expression across each Scheafer parcel
sst_expr = as.numeric(schaeff_mat[schaeff_mat$gene == 'SST',1:400])
sst_expr[is.na(sst_expr)] = 0


# transform data to prepare for parcel-wise correlations
mat_transform = as.data.frame(t(schaeff_mat[,1:400]))
colnames(mat_transform) = schaeff_mat$gene
mat_transform[is.na(mat_transform)] = 0


# check ordering
roi_part_hsq$net == gsub('7Networks',paste0(type,'_Net7'), rownames(mat_transform))

# combine gene expression with partitioned heritability
parcel_hsq_tmp = cbind(roi_part_hsq, mat_transform)

# get rid of parcels with missing gene expression
parcel_hsq_tmp = parcel_hsq_tmp[parcel_hsq_tmp$PVALB != 0,]


# remove PVALB expression outliers
upper =  mean(parcel_hsq_tmp$PVALB) + 4*sd(parcel_hsq_tmp$PVALB)
lower =  mean(parcel_hsq_tmp$PVALB) - 4*sd(parcel_hsq_tmp$PVALB)
parcel_hsq_tmp = parcel_hsq_tmp %>% filter(PVALB <= upper & PVALB >= lower)


# remove SST expression outliers
upper =  mean(parcel_hsq_tmp$SST) + 4*sd(parcel_hsq_tmp$SST)
lower =  mean(parcel_hsq_tmp$SST) - 4*sd(parcel_hsq_tmp$SST)
parcel_hsq_tmp = parcel_hsq_tmp %>% filter(SST <= upper & SST >= lower)


# PVALB heritability outliers
upper =  mean(parcel_hsq_tmp$norm_pvalb_hsq) + 4*sd(parcel_hsq_tmp$norm_pvalb_hsq)
lower =  mean(parcel_hsq_tmp$norm_pvalb_hsq) - 4*sd(parcel_hsq_tmp$norm_pvalb_hsq)
parcel_hsq_tmp = parcel_hsq_tmp %>% filter(norm_pvalb_hsq <= upper & norm_pvalb_hsq >= lower)


# SST heritability outliers
upper =  mean(parcel_hsq_tmp$norm_sst_hsq) + 4*sd(parcel_hsq_tmp$norm_sst_hsq)
lower =  mean(parcel_hsq_tmp$norm_sst_hsq) - 4*sd(parcel_hsq_tmp$norm_sst_hsq)
parcel_hsq_tmp = parcel_hsq_tmp %>% filter(norm_sst_hsq <= upper & norm_sst_hsq >= lower)

dim(parcel_hsq_tmp)

# normalized PVALB correlated to the expression of all genes
pvalb_arr    = cor(parcel_hsq_tmp$norm_pvalb_hsq, parcel_hsq_tmp[colnames(mat_transform)])
# spearman corr
pvalb_arr_sp = cor(parcel_hsq_tmp$norm_pvalb_hsq, parcel_hsq_tmp[colnames(mat_transform)], method='spearman')


# identify where PVALB falls in the distribution
#
# pearson
length(which( pvalb_arr >= pvalb_arr[colnames(pvalb_arr) == 'PVALB'] ))
print(length(which( pvalb_arr >= pvalb_arr[colnames(pvalb_arr) == 'PVALB'] ))/length(pvalb_arr))
print(cor.test(parcel_hsq_tmp$PVALB, parcel_hsq_tmp$norm_pvalb_hsq))
print(cor.test(parcel_hsq_tmp$SST, parcel_hsq_tmp$norm_pvalb_hsq))

# spearman
print(length(which( pvalb_arr_sp >= pvalb_arr_sp[colnames(pvalb_arr_sp) == 'PVALB'] ))/length(pvalb_arr_sp))
print(cor.test(parcel_hsq_tmp$PVALB, parcel_hsq_tmp$norm_pvalb_hsq, method='spearman'))


# normalized SST correlated to the expression of all genes
sst_arr = cor(parcel_hsq_tmp$norm_sst_hsq, parcel_hsq_tmp[colnames(mat_transform)])

# identify where PVALB falls in the distribution
length(which( sst_arr > sst_arr[colnames(sst_arr) == 'SST'] ))
print(length(which( sst_arr >= sst_arr[colnames(sst_arr) == 'SST'] ))/length(sst_arr))

print(cor.test(parcel_hsq_tmp$SST, parcel_hsq_tmp$norm_sst_hsq))
print(cor.test(parcel_hsq_tmp$SST, parcel_hsq_tmp$norm_sst_hsq, method='spearman'))


# control for overall heritability
out = summary(lm('scale(norm_pvalb_hsq) ~ scale(PVALB) + scale(tot_hsq)', parcel_hsq_tmp))
out


# read DFC cell type fractions
lake_dfc_dat   = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_400_7Net_expr_mat.csv'))
lake_dfc_dat_t = as.data.frame(t(lake_dfc_dat[,1:400]))
colnames(lake_dfc_dat_t) = lake_dfc_dat$gene
lake_dfc_dat_t$net = rownames(lake_dfc_dat_t)
lake_dfc_dat_t$net = gsub('7Networks','RSFA_Net7', lake_dfc_dat_t$net)


# cell type correlation to PVALB partitioned heritability
cell_hsq_df = merge(x=parcel_hsq_tmp, y=lake_dfc_dat_t, by='net')
cor(cell_hsq_df$norm_pvalb_hsq, cell_hsq_df[lake_dfc_dat$gene])
pvalb_dfc_cell = cor(cell_hsq_df$norm_pvalb_hsq, cell_hsq_df[lake_dfc_dat$gene])
pvalb_dfc_cell_df = data.frame(cor=t(pvalb_dfc_cell), cell=gsub('DFC_','',colnames(pvalb_dfc_cell)), ref='DFC')
#
cor.test(cell_hsq_df$norm_pvalb_hsq, cell_hsq_df$DFC_PVALB)
cor.test(cell_hsq_df$norm_pvalb_hsq, cell_hsq_df$DFC_PVALB, method='spearman')
cor.test(cell_hsq_df$norm_pvalb_hsq, cell_hsq_df$DFC_SST)

# SST partitioned heritability correlation
cor(cell_hsq_df$norm_sst_hsq, cell_hsq_df[lake_dfc_dat$gene])
sst_dfc_cell = cor(cell_hsq_df$norm_sst_hsq, cell_hsq_df[lake_dfc_dat$gene])
sst_dfc_cell_df = data.frame(cor=t(sst_dfc_cell), cell=gsub('DFC_','',colnames(sst_dfc_cell)), ref='DFC')



# read VIS cell type fractions
lake_vis_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_VIS_',nparcel,'_',net_num,'Net_expr_mat.csv'))
lake_vis_dat_t = as.data.frame(t(lake_vis_dat[,1:400]))
colnames(lake_vis_dat_t) = lake_vis_dat$gene
lake_vis_dat_t$net = rownames(lake_vis_dat_t)
lake_vis_dat_t$net = gsub('7Networks','RSFA_Net7', lake_vis_dat_t$net)

# cell type correlation
vis_cell_hsq_df = merge(x=parcel_hsq_tmp, y=lake_vis_dat_t, by='net')
cor(vis_cell_hsq_df$norm_pvalb_hsq, vis_cell_hsq_df[lake_vis_dat$gene])
pvalb_vis_cell = cor(vis_cell_hsq_df$norm_pvalb_hsq, vis_cell_hsq_df[lake_vis_dat$gene])
pvalb_vis_cell_df = data.frame(cor=t(pvalb_vis_cell), cell=gsub('VIS_','',colnames(pvalb_vis_cell)), ref='VIS')

cor(vis_cell_hsq_df$norm_sst_hsq, vis_cell_hsq_df[lake_vis_dat$gene])
sst_vis_cell    = cor(vis_cell_hsq_df$norm_sst_hsq, vis_cell_hsq_df[lake_vis_dat$gene])
sst_vis_cell_df = data.frame(cor=t(sst_vis_cell), cell=gsub('VIS_','',colnames(sst_vis_cell)), ref='VIS')
#
cor.test(vis_cell_hsq_df$norm_pvalb_hsq, vis_cell_hsq_df$VIS_PVALB)
cor.test(vis_cell_hsq_df$norm_pvalb_hsq, vis_cell_hsq_df$VIS_PVALB, method='spearman')
cor.test(vis_cell_hsq_df$norm_pvalb_hsq, vis_cell_hsq_df$VIS_SST)


# VIS/DFC cells correlated to PVALB partitioned heritability
pvalb_cell_df = rbind(pvalb_vis_cell_df, pvalb_dfc_cell_df)

# VIS/DFC cells correlated to SST partitioned heritability
sst_cell_df   = rbind(sst_vis_cell_df, sst_dfc_cell_df)


sst_cell_df = rbind(sst_cell_df, data.frame(cor=NA, cell='In2', ref='DFC'))
sst_cell_df = rbind(sst_cell_df, data.frame(cor=NA, cell='In3', ref='VIS'))
sst_cell_df = rbind(sst_cell_df, data.frame(cor=NA, cell='Ex2', ref='VIS'))


# determine cell type ordering for the plot
avg_cor = sst_cell_df %>% group_by(cell) %>% summarise(x=mean(cor, na.rm=T))
sst_cell_df$cell = factor(sst_cell_df$cell, levels=avg_cor$cell[order(avg_cor$x)])
pvalb_cell_df$ref = factor(pvalb_cell_df$ref, levels=c('DFC','VIS'))


figure_out = paste0(base_dir,'/figures/PaperFig_celltype_CELLTYPES_to_SST_partitioned_heritability.pdf')
figure_out
pdf(figure_out, width=2.5, height=1.5)
ggplot(data=sst_cell_df, aes(x=cell, y=cor, fill=ref)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits=c(-.4,.4), breaks = seq(-.4, .4, by = .2)) +
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


# VIS/DFC cells correlated to SST partitioned heritability
pvalb_cell_df = rbind(pvalb_cell_df, data.frame(cor=NA, cell='In2', ref='DFC'))
pvalb_cell_df = rbind(pvalb_cell_df, data.frame(cor=NA, cell='In3', ref='VIS'))
pvalb_cell_df = rbind(pvalb_cell_df, data.frame(cor=NA, cell='Ex2', ref='VIS'))


# determine cell type ordering for the plot
avg_cor = pvalb_cell_df %>% group_by(cell) %>% summarise(x=mean(cor, na.rm=T))
pvalb_cell_df$cell = factor(pvalb_cell_df$cell, levels=avg_cor$cell[order(avg_cor$x)])
pvalb_cell_df$ref = factor(pvalb_cell_df$ref, levels=c('DFC','VIS'))

figure_out = paste0(base_dir,'/figures/PaperFig_celltype_CELLTYPES_to_PVALB_partitioned_heritability.pdf')
figure_out
pdf(figure_out, width=2.5, height=1.5)
ggplot(data=pvalb_cell_df, aes(x=cell, y=cor, fill=ref)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
    #geom_bar(stat="identity", color="black", position=position_dodge(), bar_spacing=1) +
    #geom_errorbar(aes(ymin=norm_expr-expr_se, ymax=norm_expr+expr_se), width=.4, position=position_dodge(.9)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits=c(-.4,.4), breaks = seq(-.4, .4, by = .2)) +
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









# Plot Partioned Heritability versus gene expression
# ------------
# network colors
rgb_arr = c(rgb(42/255,204/255,160/255),
             rgb(124/255,18/255,133/255),
             rgb(230/255,147/255,36/255),
             rgb(70/255,130/255,182/255),
             rgb(205/255,61/255,84/255),
             rgb(219/255,248/255,167/255),
             rgb(122/255,135/255,50/255))

cor.test(parcel_hsq_tmp$norm_pvalb_hsq, parcel_hsq_tmp$PVALB)
cor.test(parcel_hsq_tmp$norm_sst_hsq, parcel_hsq_tmp$SST)

parcels = read.csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_n9713_cluster_7.csv'))
parcel_hsq_plot = merge(x=parcel_hsq_tmp, by.x='net', y=parcels, by.y='parcels')

# plot expression by heritability
out_file = paste0(base_dir, '/figures/PaperFig_clust7_schaeff400_pvalbHSQ_parcelwise.pdf')
CairoPDF(width=2, height=2, file=out_file)
p = ggplot(parcel_hsq_plot, aes(PVALB, norm_pvalb_hsq)) +
  geom_point(aes(colour = factor(clusters)), size = 1.5, alpha=1) +
  scale_color_manual(values = rgb_arr) +
  theme_classic() +
  xlab('PVALB expression') +
  ylab(paste0('h2snp:pvalb / h2snp:total')) +
  expand_limits(x = c(-3, 3)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-3, 3, by = 1)) +
  expand_limits(y = c(0, .5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, .5, by = .1)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), legend.position = "none")
print(p)
dev.off()


# plot expression by heritability
out_file = paste0(base_dir, '/figures/PaperFig_clust7_schaeff400_sstHSQ_parcelwise.pdf')
CairoPDF(width=2, height=2, file=out_file)
p = ggplot(parcel_hsq_plot, aes(SST, norm_sst_hsq)) +
  geom_point(aes(colour = factor(clusters)), size = 1.5, alpha=1) +
  scale_color_manual(values = rgb_arr) +
  theme_classic() +
  xlab('PVALB expression') +
  ylab(paste0('h2snp:sst / h2snp:total')) +
  expand_limits(x = c(-3, 3)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-3, 3, by = 1)) +
  expand_limits(y = c(0, .5)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, .5, by = .1)) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), legend.position = "none")
print(p)
dev.off()





# Density plot of PVALB-to-PVALB Heritability correlation, relative to all genes
cor_df    = data.frame(cors=as.numeric(pvalb_arr), gene=colnames(pvalb_arr))
pvalb_cor = cor_df$cors[cor_df$gene == 'PVALB']
pop_mean  = mean(cor_df$cors[cor_df$gene != 'PVALB'])
pop_sd    = sd(cor_df$cors[cor_df$gene != 'PVALB'])

length(which(cor_df$cors >= pvalb_cor))
length(cor_df$cors)


pdf('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_rsfa_pvalb_h2_to_pvalbexpr_distribution.pdf',width=5,height=3)
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
    scale_x_continuous(expand = c(0, 0), limits=c(-.5,.5), breaks = seq(-.5,.5,.25)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = pvalb_cor, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()




# Density plot of SST-to-SST Heritability correlation, relative to all genes
cor_df   = data.frame(cors=as.numeric(sst_arr), gene=colnames(sst_arr))
sst_cor  = cor_df$cors[cor_df$gene == 'SST']
pop_mean = mean(cor_df$cors[cor_df$gene != 'SST'])
pop_sd   = sd(cor_df$cors[cor_df$gene != 'SST'])

length(which(cor_df$cors <= sst_cor))
length(cor_df$cors)


pdf('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_rsfa_sst_h2_to_sstexpr_distribution.pdf',width=5,height=3)
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
    scale_x_continuous(expand = c(0, 0), limits=c(-.3,.3), breaks = round(seq(-.3,.3,.1),2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = sst_cor, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()





parcel_hsq = cur_dat
nparcel = '400'
net_num = '7'
plot_matlab(values=as.numeric(parcel_hsq$norm_pvalb_hsq), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_rsfa_pvalb_norm_hsq.dscalar.nii'), parcel_num=nparcel, net_num=net_num)


# Figure 5 - Plot PVALB Partitioned heritability
# --------
map2color = function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


# get rid of extreme hs2part values
parcel_hsq$norm_pvalb_hsq[parcel_hsq$norm_pvalb_hsq > .75] <- .000001 # so they have the most negative color in wb_view, rather than being blank

# create custom pvalb pallete
pvalb_pal = colorRampPalette(c("white", '#EBF0B5', '#46B6C2', '#283B8E'))(20)
# color scale for figure 1
setEPS()
postscript(width=2, height=4,file=paste0(base_dir, '/figures/colorbar_pvalb.eps'))
image(x,y,z,col=parcel_hsq$norm_pvalb_hsq, axes=FALSE,xlab="",ylab="")
dev.off()

label_file = paste0(base_dir, '/figures/label_file.txt')
sink(label_file)
color_labels = map2color(parcel_hsq$norm_pvalb_hsq, pvalb_pal, limits=c(0,.25))
label_array = rep(0,400)
col_ct = 1
for (col in unique(color_labels)){
  label_array[color_labels == col] <- col_ct
  cat(as.character(col_ct))
  cat('\n')
  cat(paste(as.character(col_ct), as.character(col2rgb(col)[1]), as.character(col2rgb(col)[2]), as.character(col2rgb(col)[3]), 255))
  cat('\n')
  col_ct <- col_ct + 1
}
sink()

save_path = paste0(base_dir, '/figures/PaperFig_clust7_noscale_ahba_pvalb_partial_hsq_rsfa_0.25.dscalar.nii')
plot_matlab(values=parcel_hsq$norm_pvalb_hsq, out_path=save_path, parcel_num=nparcel, net_num=net_num)

save_path = paste0(base_dir, '/figures/PaperFig_clust7_noscale_ahba_pvalb_partial_hsq_rsfa_0.25.dlabel.nii')
plot_matlab(values=label_array, out_path=save_path, parcel_num=nparcel, net_num=net_num)
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_file, ' ', save_path))




# Figure 5 - Plot SST Partitioned heritability
# --------
sst_pal <- colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D'))(20)
setEPS()

postscript(width=2, height=4,file=paste0(base_dir, '/figures/colorbar_sst.eps'))
image(x,y,z,col=parcel_hsq$norm_sst_hsq,axes=FALSE,xlab="",ylab="")
dev.off()

parcel_hsq$norm_sst_hsq[parcel_hsq$norm_sst_hsq > .90] <- .000001 # so they have the most negative color in wb_view, rather than being blank

label_file   = paste0(base_dir, '/figures/label_file.txt')
sink(label_file)
color_labels = map2color(parcel_hsq$norm_sst_hsq, sst_pal, limits=c(0,.25))
label_array  = rep(0,400)
col_ct = 1
for (col in unique(color_labels)){
  label_array[color_labels == col] <- col_ct
  cat(as.character(col_ct))
  cat('\n')
  cat(paste(as.character(col_ct), as.character(col2rgb(col)[1]), as.character(col2rgb(col)[2]), as.character(col2rgb(col)[3]), 255))
  cat('\n')
  col_ct <- col_ct + 1
}
sink()

save_path = paste0(base_dir, '/figures/clust7_ahba_sst_partial_hsq_rsfa_0.25.dlabel.nii')
plot_matlab(values=label_array, out_path=save_path, parcel_num=nparcel, net_num=net_num)
system(paste0('/nexsan/apps/hpc/Apps/ConnectomeWorkbench/1.2.3/workbench/bin_rh_linux64/wb_command -cifti-label-import ', save_path, ' ', label_file, ' ', save_path))







