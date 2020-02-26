
library(tidyverse)
library(Cairo)

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
base_dir      = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
rsfa_gwas_dir = paste0(base_dir, '/data/ukb/imputed/rsfa_gwas')
ldsc_part_dir = paste0(base_dir, '/data/ukb/imputed/ldsc_part')


# read schaeffer averaged expression data
schaeff_mat  = read_csv(paste0(base_dir, '/data/ahba/schaeffer400_7Net_expr_mat.csv'))
expr_parcels = colnames(schaeff_mat)[1:400]


# prepare LDSC partitioned h2 values for plotting
nparcel = '400'
net_num = '7'
schaeffer_parcels = read.table('/gpfs/milgram/project/holmes/HOLMES_UKB/atlas/Schaefer2018_LocalGlobal/Parcellations/MNI/Schaefer2018_400Parcels_7Networks_order.txt')
tmp_plot_df = data.frame(idx=1:400, parcel_name=gsub('7Networks','RSFA_Net7', schaeffer_parcels$V2))


# transform data to prepare for parcel-wise correlations
mat_transform = as.data.frame(t(schaeff_mat[,1:400]))
colnames(mat_transform) = schaeff_mat$gene
mat_transform[is.na(mat_transform)] = 0

# change label names for merge
mat_transform$network = gsub('7Networks','RSFA_Net7', rownames(mat_transform))


# schaefer parcel names (important for mapping)
analysis_dir = paste0(base_dir, '/data/ukb/rsfa/pheno_RSFA')
parcel_names = read_csv(paste0(analysis_dir, '/RSFA_brit_n9713_400parcel_colnames.txt'), col_names=FALSE)


# Stage 1: LDSC OVERALL
#
# read the LDSC heritability estimates
parcel_h2_df = NULL
for (parcel in 1:400){
    write(parcel,'')
    parcel_str = str_pad(parcel,3,'left','0')

    in_file = paste0(rsfa_gwas_dir, '/rsfa_brit_n9713_parcel_',parcel_str,'_400_fastGWA_eur_w_ld_chr_ldsc_h2.log')
    con     = file(in_file, 'r')
    while ( TRUE ){
        line = readLines(con,n=1)
        if (length(line) == 0){
            break
        }
        if (length(grep('Total Observed',line)) == 1){
            dat_1 = strsplit(line, 'h2:')[[1]][[2]]
            h2    = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[2]]
            h2_se = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[3]]
        }
    }
    out_row = data.frame(parcel=parcel, parcel_name=parcel_names$X1[parcel], h2=as.numeric(h2), h2_se=as.numeric(h2_se))
    parcel_h2_df = rbind(parcel_h2_df, out_row)
}


# read the GCTA based heritability estimatse
gcta_h2      = read_csv(paste0(base_dir, '/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv'))
gcta_ldsc_h2 = merge(x=gcta_h2, y=parcel_h2_df, by.x='net', by.y='parcel_name')
cor.test(gcta_ldsc_h2$tot_hsq, gcta_ldsc_h2$h2, method='spearman')

ldsc_gcta_plot = ggplot(data=gcta_ldsc_h2, aes(x=tot_hsq, y=h2)) +
                        geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = .5) +
                        theme_classic() +
                        geom_smooth(method = lm, se = FALSE, color='black', linetype='dashed', fullrange=T) +
                        scale_x_continuous(limits=c(0,.4), breaks=seq(0,.4,.1), expand=c(0,0)) +
                        scale_y_continuous(limits=c(0,.4), breaks=seq(0,.4,.1), expand=c(0,0)) +
                        xlab('GCTA Total Heritability') +
                        ylab('LDSC Total Heritability')
fig_out = paste0(base_dir, '/figures/gcta_ldsc_total_heritability_corr_plot.pdf')

CairoPDF(fig_out, width=4, height=4)
print(ldsc_gcta_plot)
dev.off()
fig_out


# plot total heritability
total_plot_df  = merge(x=tmp_plot_df, y=parcel_h2_df, by.x='parcel_name', by.y='parcel_name')
total_plot_df  = total_plot_df[order(total_plot_df$idx),]
# wb_view has trouble plotting zeros, so just make them very tiny
total_plot_df$h2[total_plot_df$h2 < 0] = 0.0001
out_path = paste0(base_dir, '/figures/surface_plots/LDSC_PaperFig_total_heritability.dscalar.nii')
plot_matlab(values= total_plot_df$h2, out_path=out_path, parcel_num=nparcel, net_num=net_num)



# Stage 2: LDSC PARTITIONED (PVALB/SST)
#
# Read partitioned heritability estimates from LDSC
part_parcel_h2_df = NULL
# read the LDSC heritability estimates
for (gene in c('SST','PVALB')){
    for (parcel in 1:400){
        write(parcel,'')
        parcel_str = str_pad(parcel,3,'left','0')

        in_file        = paste0(ldsc_part_dir, '/GCTA_',gene,'_ahba_top500_ensembl_rsfa_brit_n9713_parcel_',parcel_str,'_400.results')
        part_herit_in  = read.table(in_file, header=T)
        out_row        = part_herit_in[part_herit_in$Category == 'L2_0',]
        out_row$parcel = parcel
        out_row$net    = parcel_names$X1[parcel]
        out_row$gene   = gene
        part_parcel_h2_df = rbind(part_parcel_h2_df, out_row)
    }
}
pvalb_parcel_h2_df = part_parcel_h2_df %>% filter(gene== 'PVALB')
sst_parcel_h2_df   = part_parcel_h2_df %>% filter(gene== 'SST')

# unlike GCTA, LDSC doesnt constrain h2snp to be above zero
pvalb_parcel_h2_df$Prop._h2[pvalb_parcel_h2_df$Prop._h2 < 0] = 0
sst_parcel_h2_df$Prop._h2[sst_parcel_h2_df$Prop._h2 < 0] = 0


# merge LDSC parcel-wise partitioned heritability with AHBA expression data
pvalb_h2_w_expr_df = merge(x=mat_transform, by.x='network', y=pvalb_parcel_h2_df, by.y='net')
sst_h2_w_expr_df   = merge(x=mat_transform, by.x='network', y=sst_parcel_h2_df, by.y='net')


# get rid of parcels with no expression dat
pvalb_h2_w_expr_df = pvalb_h2_w_expr_df %>% filter(pvalb_h2_w_expr_df$PVALB != 0)
sst_h2_w_expr_df   = sst_h2_w_expr_df %>% filter(sst_h2_w_expr_df$SST != 0)


# no SST/PVALB expression outliers
which(abs(scale(pvalb_h2_w_expr_df$PVALB)) > 4)
which(abs(scale(sst_h2_w_expr_df$SST)) > 4)

# some heritability outliers
which(abs(scale(pvalb_h2_w_expr_df$Prop._h2)) > 4)
which(abs(scale(sst_h2_w_expr_df$Prop._h2)) > 4)

# get rid of partitioned heritability outliers
pvalb_h2_w_expr_df = pvalb_h2_w_expr_df[which(abs(scale(pvalb_h2_w_expr_df$Prop._h2)) < 4),]
sst_h2_w_expr_df   = sst_h2_w_expr_df[which(abs(scale(sst_h2_w_expr_df$Prop._h2)) < 4),]

# data frame with gene-wise correlation to PVALB heritability map
all_pvalb_gene_cors = cor(pvalb_h2_w_expr_df[schaeff_mat$gene], pvalb_h2_w_expr_df$Prop._h2, method='spearman')
pvalb_ldsc_cor_df   = data.frame(cor=all_pvalb_gene_cors, gene=rownames(all_pvalb_gene_cors))
pvalb_ldsc_cor_df   = pvalb_ldsc_cor_df[rev(order(pvalb_ldsc_cor_df$cor)),]

pvalb_ldsc_cor_df[which(pvalb_ldsc_cor_df$gene == 'PVALB'),]
which(pvalb_ldsc_cor_df$gene == 'PVALB')

# data frame with gene-wise correlation to SST heritability map
all_sst_gene_cors = cor(sst_h2_w_expr_df[schaeff_mat$gene], sst_h2_w_expr_df$Enrichment, method='spearman')
sst_ldsc_cor_df   = data.frame(cor=all_sst_gene_cors, gene=rownames(all_sst_gene_cors))
sst_ldsc_cor_df   = sst_ldsc_cor_df[rev(order(sst_ldsc_cor_df$cor)),]

sst_ldsc_cor_df[which(sst_ldsc_cor_df$gene == 'SST'),]
which(sst_ldsc_cor_df$gene == 'SST')


# PVALB - plot values on cortical surface
pvalb_plot_df  = merge(x=tmp_plot_df, y=pvalb_parcel_h2_df, by.x='parcel_name', by.y='net')
pvalb_plot_df  = pvalb_plot_df[order(pvalb_plot_df$idx),]
out_path = paste0(base_dir, '/figures/surface_plots/Phase1_LDSC_PaperFig_rsfa_pvalb_prop_h2.dscalar.nii')
plot_matlab(values= pvalb_plot_df$Enrichment, out_path=out_path, parcel_num=nparcel, net_num=net_num)

# SST - plot values on cortical surface
plot_df  = merge(x=tmp_plot_df, y=sst_parcel_h2_df, by.x='parcel_name', by.y='net')
plot_df  = plot_df[order(plot_df$idx),]
out_path = paste0(base_dir, '/figures/surface_plots/Phase1_LDSC_PaperFig_rsfa_sst_prop_h2.dscalar.nii')
plot_matlab(values= plot_df$Enrichment, out_path=out_path, parcel_num=nparcel, net_num=net_num)




# read the GCTA based heritability estimatse
gcta_h2      = read_csv(paste0(base_dir, '/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv'))
head(gcta_h2)
gcta_ldsc_h2 = merge(x=gcta_h2, y=pvalb_parcel_h2_df, by.x='net', by.y='net')

# get rid of any heritability outliers
gcta_ldsc_h2 = gcta_ldsc_h2[which(abs(scale(gcta_ldsc_h2$norm_pvalb_hsq)) < 4),]
gcta_ldsc_h2 = gcta_ldsc_h2[which(abs(scale(gcta_ldsc_h2$Prop._h2)) < 4),]

cor.test(gcta_ldsc_h2$norm_pvalb_hsq, gcta_ldsc_h2$Enrichment, method='spearman')

ldsc_gcta_plot = ggplot(data=gcta_ldsc_h2, aes(x=norm_pvalb_hsq, y=Prop._h2)) +
                        geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = .5) +
                        theme_classic() +
                        geom_smooth(method = lm, se = FALSE, color='black', linetype='dashed', fullrange=T) +
                        scale_x_continuous(limits=c(0,.35), breaks=seq(0,.35,.05), expand=c(0,0)) +
                        scale_y_continuous(limits=c(0,.35), breaks=seq(0,.35,.05), expand=c(0,0)) +
                        xlab('GCTA PVALB partitioned heritability') +
                        ylab('LDSC PVALB partitioned heritability')
fig_out = paste0(base_dir, '/figures/gcta_ldsc_PVALB_corr_plot.pdf')

CairoPDF(fig_out, width=4, height=4)
print(ldsc_gcta_plot)
dev.off()
fig_out


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

cor.test(pvalb_h2_w_expr_df$Prop._h2, pvalb_h2_w_expr_df$PVALB, method='spearman')

parcels = read.csv(paste0(base_dir, '/data/ukb/rsfa/RSFA_n9713_cluster_7.csv'))
parcel_hsq_plot = merge(x=pvalb_h2_w_expr_df, by.x='network', y=parcels, by.y='parcels')

# plot expression by heritability
out_file = paste0(base_dir, '/figures/PaperFig_clust7_schaeff400_LDSC_pvalbHSQ_parcelwise.pdf')
CairoPDF(width=3, height=3, file=out_file)
p = ggplot(parcel_hsq_plot, aes(x=PVALB, y=Prop._h2)) +
              geom_point(aes(colour = factor(clusters)), size = 1.5, alpha=1) +
              scale_color_manual(values = rgb_arr) +
              theme_classic() +
              xlab('PVALB expression') +
              ylab(paste0('h2snp:pvalb / h2snp:total')) +
              expand_limits(x = c(-3, 3)) +
              scale_x_continuous(expand = c(0, 0), breaks = seq(-3, 3, by = 1)) +
              expand_limits(y = c(0, .5)) +
              scale_y_continuous(expand = c(0, 0), limits=c(0,.4), breaks = seq(0, .4, by = .2)) +
              geom_smooth(method = "lm", se = FALSE, linetype='dashed', color='black', fullrange=T) +
              theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), legend.position = "none")
print(p)
dev.off()




# SST partitioned heritability to SST expression
cor.test(sst_h2_w_expr_df$Prop._h2, sst_h2_w_expr_df$SST, method='spearman')
sst_parcel_hsq_plot = merge(x=sst_h2_w_expr_df, by.x='network', y=parcels, by.y='parcels')

# plot expression by heritability
out_file = paste0(base_dir, '/figures/PaperFig_clust7_schaeff400_LDSC_sstHSQ_parcelwise.pdf')
CairoPDF(width=3, height=3, file=out_file)
p = ggplot(sst_parcel_hsq_plot, aes(x=SST, y=Prop._h2)) +
              geom_point(aes(colour = factor(clusters)), size = 1.5, alpha=1) +
              scale_color_manual(values = rgb_arr) +
              theme_classic() +
              xlab('SST expression') +
              ylab(paste0('h2snp:pvalb / h2snp:total')) +
              expand_limits(x = c(-3, 3)) +
              scale_x_continuous(expand = c(0, 0), breaks = seq(-3, 3, by = 1)) +
              expand_limits(y = c(0, .5)) +
              scale_y_continuous(expand = c(0, 0), limits=c(0,.4), breaks = seq(0, .4, by = .2)) +
              geom_smooth(method = "lm", se = FALSE, linetype='dashed', color='black', fullrange=T) +
              theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"), legend.position = "none")
print(p)
dev.off()




# Stage 3: LDSC 7 CLUSTER OVERALL
#
# read the LDSC heritability estimates
clust_names   = read.csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ukb/rsfa/pheno_RSFA/RSFA_brit_n9713_cluster7_colnames.txt')
clust_names$x = as.character(clust_names$x)
cluster_h2_df = NULL
for (parcel in 1:7){
    write(parcel,'')
    parcel_str = str_pad(parcel,3,'left','0')
    in_file    = paste0(rsfa_gwas_dir, '/rsfa_brit_n9713_cluster_',parcel_str,'_7_fastGWA_eur_w_ld_chr_ldsc_h2.log')

    con = file(in_file, 'r')
    while ( TRUE ){
        line = readLines(con,n=1)
        if (length(line) == 0){
            break
        }
        if (length(grep('Total Observed',line)) == 1){
            dat_1 = strsplit(line, 'h2:')[[1]][[2]]
            h2    = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[2]]
            h2_se = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[3]]
        }
    }
    out_row = data.frame(parcel=parcel, parcel_name=clust_names$x[parcel], h2=as.numeric(h2), h2_se=as.numeric(h2_se))
    cluster_h2_df = rbind(cluster_h2_df, out_row)
}


# order the networks for plotting
positions = c('limbicA','limbicB','cingulo-opercular','prefrontal','parietal','motor','visual')
cluster_h2_df$region = factor(cluster_h2_df$parcel_name, levels=positions)
color_arr = c('#DBF8A7','#7A8732','#2ACCA0','#CD3D54','#E69324','#4682B6','#7C1285')


# Figure 5b replication
fig_out = paste0(base_dir, '/figures/PaperFig_LDSC_clust7_noscale_hsq.pdf')
fig_out
CairoPDF(fig_out, width=5, height=3)
p1 = ggplot(cluster_h2_df, aes(x=region, y=h2)) +
                geom_bar(aes(fill = region), stat="identity", color="black", position=position_dodge()) +
                geom_errorbar(aes(fill = region, ymin=h2-h2_se, ymax=h2+h2_se), width=.4, position=position_dodge(.9)) +
                theme_classic() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                scale_y_continuous(breaks = seq(0,.35, .05), expand = c(0, 0), limits=c(0,.35)) +
                theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black")) +
      scale_fill_manual(values = color_arr)

print(p1)
dev.off()



# Stage 4: LDSC 7 CLUSTER PARTITIONED (SST/PVALB)
#
clust_part_h2_df = NULL
for (gene in c('PVALB','SST')){
    for (parcel in 1:7){
        write(parcel,'')
        parcel_str = str_pad(parcel,3,'left','0')

        in_file        = paste0(ldsc_part_dir, '/GCTA_',gene,'_ahba_top500_ensembl_rsfa_brit_n9713_cluster_',parcel_str,'_7.results')
        part_herit_in  = read.table(in_file, header=T)
        out_row        = part_herit_in[part_herit_in$Category == 'L2_0',]
        out_row$parcel = parcel
        out_row$net    = clust_names$x[parcel]
        out_row$gene   = gene
        clust_part_h2_df = rbind(clust_part_h2_df, out_row)
    }
}
clust_part_h2_df$Prop._h2[clust_part_h2_df$Prop._h2 < 0] = 0

ggplot(clust_part_h2_df, aes(x=net, y=Enrichment, fill=gene)) +
    geom_bar(stat='identity', position='dodge')


# FIgure 5b
clust_part_h2_df$net = factor(clust_part_h2_df$net, levels=c('limbicA','limbicB','cingulo-opercular','parietal','prefrontal','motor','visual'))
fig_out = paste0(base_dir, '/figures/PaperFig_LDSC_clust7_noscale_partitioned_hsq_interneuron.pdf')
fig_out
pdf(fig_out, width=3, height=3)
p1 = ggplot(clust_part_h2_df, aes(x=net, y=Prop._h2)) +
              geom_bar(aes(fill = gene), stat="identity", color="black", position=position_dodge()) +
              geom_errorbar(aes(fill = gene, ymin=Prop._h2-Prop._h2_std_error, ymax=Prop._h2+Prop._h2_std_error), width=.4,  position=position_dodge(.9)) +
              theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
              scale_y_continuous(breaks=seq(-0.1,.25,.05), expand = c(0, 0), limits=c(-0.1,.25)) +
              theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
print(p1)
dev.off()





# read the LDSC heritability estimates
parcel_h2_df = NULL
for (parcel in 1:400){
    write(parcel,'')
    parcel_str = str_pad(parcel,3,'left','0')

    in_file = paste0(rsfa_gwas_dir, '/rsfa_brit_n9713_parcel_',parcel_str,'_400_fastGWA_eur_w_ld_chr_ldsc_h2.log')
    con     = file(in_file, 'r')
    while ( TRUE ){
        line = readLines(con,n=1)
        if (length(line) == 0){
            break
        }
        if (length(grep('Total Observed',line)) == 1){
            dat_1 = strsplit(line, 'h2:')[[1]][[2]]
            h2    = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[2]]
            h2_se = strsplit(gsub('[(]|[)]','',dat_1), ' ')[[1]][[3]]
        }
    }
    out_row = data.frame(parcel=parcel, parcel_name=parcel_names$X1[parcel], h2=as.numeric(h2), h2_se=as.numeric(h2_se))
    parcel_h2_df = rbind(parcel_h2_df, out_row)
}


# read the GCTA based heritability estimatse
gcta_h2      = read_csv(paste0(base_dir, '/data/ukb/rsfa/PaperData_all_RSFA_maf05_partition_parcel_hsq.csv'))
gcta_ldsc_h2 = merge(x=gcta_h2, y=parcel_h2_df, by.x='net', by.y='parcel_name')
cor.test(gcta_ldsc_h2$tot_hsq, gcta_ldsc_h2$h2)
gcta_ldsc_h2 = gcta_ldsc_h2[gcta_ldsc_h2$h2 > 0,]


# plot correspondence between the two methods
h2_plot = ggplot(gcta_ldsc_h2, aes(x=h2, y=tot_hsq), alpha=0.8) +
                    geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 1) +
                    theme_classic() +
                    scale_x_continuous(limits=c(0,.36), expand=c(0,0), breaks=seq(0,.36, .12)) +
                    scale_y_continuous(limits=c(0,.36), expand=c(0,0), breaks=seq(0,.36, .12)) +
                    geom_smooth(method='lm', formula= y~x, color='black', fullrange=T, se=F, linetype='dashed') +
                    theme(text = element_text(size=12, color='black')) +
                    xlab('LDSC heritability') +
                    ylab('GCTA heritability')

out_path = paste0(base_dir, '/figures/gcta_ldsc_correspondence.pdf')
out_path
CairoPDF(out_path, width=3, height=3)
print(h2_plot)
dev.off()





