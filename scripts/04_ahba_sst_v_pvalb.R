library(R.matlab)
library(tidyverse)
library(stats)
library(Cairo)
source('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/scripts/R_function_library.R')


# load sample information
base_dir   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
figure_dir = paste0(base_dir, '/figures')


# read previously processed AHBA data
load_me = paste0(base_dir, '/data/ahba/donorDat_obj.Rdata')
load_me
load(file=load_me, verbose=T)


# sample information file with vertex mapping
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))


# remove samples greater than 4mm away from nearest vertex
usable_ctx_samples = sample_info[which(abs(sample_info$mm_to_surf) < 4),]


donor_nums = c('9861', '10021', '12876', '14380', '15496', '15697')


# z-transform cortical data separately for each donor
ctx_data_scale = matrix()
ctx_samp_all   = NULL
for (donor in donor_nums){
    print(paste0('Mean and variance normalizing cortical data for : ', donor))
    
    # subset to current donor, cortex only data
    ctx_samples = usable_ctx_samples %>% filter(brain == donor) %>% filter(top_level == 'CTX')
    ctx_idxs    = which(donorDat$samples$well_id %in% ctx_samples$well_id)
    ctx_micro   = donorDat$micro[,ctx_idxs]
    ctx_samp    = donorDat$samples[ctx_idxs,]

    # within subject scaling
    ctx_micro_scale = as.data.frame(t(apply(ctx_micro, 1, scale, center=TRUE, scale=TRUE)))
    rownames(ctx_micro_scale) = donorDat$probes$gene_symbol

    # build overall matrix with data from all subjects
    ctx_data_scale = cbind(ctx_data_scale, ctx_micro_scale)
    ctx_samp_all   = rbind(ctx_samp_all, ctx_samp)
}
ctx_data_scale$ctx_data_scale = NULL
ctx_samp$ctx_samp = NULL

# save scales cortical AHBA expression data
save(ctx_data_scale, ctx_samp_all, file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'))

# load from previous save-point
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)



# correlate cortical SST/PVALB
donor      = '9861'
sst_idx    = which(donorDat$probes$gene_symbol == 'SST')
pvalb_idx  = which(donorDat$probes$gene_symbol == 'PVALB')
sst_expr   = as.numeric(ctx_data_scale[rownames(ctx_data_scale) == 'SST',])
pvalb_expr = as.numeric(ctx_data_scale[rownames(ctx_data_scale) == 'PVALB',])

# pearson
sst_pvalb_cor = cor(sst_expr, pvalb_expr)

# spearman (sst/pvalb expression are fairly normally distributed, but good to show consistency with non-parametric stats)
sst_pvalb_spearman_cor = cor(sst_expr, pvalb_expr, method='spearman')


# cortical correlations of every gene to sst/pvalb
sst_cor      = t(cor(sst_expr, t(ctx_data_scale)))
pvalb_cor    = t(cor(pvalb_expr, t(ctx_data_scale)))

# same, but with spearman correlations
sst_sp_cor   = t(cor(sst_expr, t(ctx_data_scale), method='spearman'))
pvalb_sp_cor = t(cor(pvalb_expr, t(ctx_data_scale), method='spearman'))


# SST - spearman
length(which(sst_sp_cor <= sst_pvalb_spearman_cor))
length(which(sst_sp_cor <= sst_pvalb_spearman_cor))/length(sst_sp_cor)

# SST - pearson
length(which(sst_cor <= sst_pvalb_cor))
length(which(sst_cor <= sst_pvalb_cor))/length(pvalb_cor)

# PVALB - spearman
length(which(pvalb_sp_cor <= sst_pvalb_spearman_cor))
length(which(pvalb_sp_cor <= sst_pvalb_spearman_cor))/length(pvalb_sp_cor)

# PVALB - pearson
length(which(pvalb_cor <= sst_pvalb_cor))
length(which(pvalb_cor <= sst_pvalb_cor))/length(pvalb_cor)


# write pvalb gene list
sst_df      = data.frame(gene_symbol=rownames(sst_cor), sst_cor=as.numeric(sst_cor))
sst_write   = merge(x=donorDat$probes, y=sst_df, by='gene_symbol')
pvalb_df    = data.frame(gene_symbol=rownames(pvalb_cor), pvalb_cor=as.numeric(pvalb_cor))
pvalb_write = merge(x=donorDat$probes, y=pvalb_df, by='gene_symbol')
#
write_csv(pvalb_write, paste0(base_dir, '/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
write_csv(sst_write, paste0(base_dir, '/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))



# compare the cortical PVALB/SST correlation to all other possible pairwise correlations
all_cors     = cor(t(ctx_data_scale))
# get unique correlations in matrix
all_cor_vals = all_cors[upper.tri(all_cors)]

# number of correlations less than the SST/PVALB one
length(which(all_cor_vals <= sst_pvalb_cor))/length(all_cor_vals)
df = data.frame(cors=all_cor_vals)

# figure shows where the SST/PVALB correlation falls relative to all other pairwise comparisons
CairoPDF('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_SSTvPVALB_relative_to_all_other_pairwise.pdf',family='Arial',width=5,height=3)
ggplot(df, aes(cors)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('All pairwise cortical spatial correlations')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=12),
          axis.text.y = element_text(colour="black", size=12),
          axis.title.x = element_text(colour="black", size=12),
          axis.title.y = element_text(colour="black", size=12)) +
    ylab("Density") +
    expand_limits(y=c(0,2)) +
    scale_x_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1,1,.5)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = sst_pvalb_cor, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()



# SST to PVALB correlation for each subcortical region
ahba_probes = donorDat$probes

# dictionary with region info
top_level        = read_csv(paste0(base_dir, '/ref_files/top_level_categories.csv'))
reg_df           = NULL
reg_df['CTX']    = 'Cortex'
reg_df['THAL']   = 'Thalamus'
reg_df['HIPP']   = 'Hippocampus'
reg_df['STR']    = 'Striatum'
reg_df['VTA_SN'] = 'VTA & SN'
reg_df['AMYG']   = 'Amygdala'
reg_df['HYPO']   = 'Hypothalamus'
reg_df['GP']     = 'Globus Pallidus'

# (pre-allocate) matrix to store correlation of every gene to SST
cor_table_SST           = data.frame(matrix(nrow=nrow(ahba_probes), ncol=length(reg_df), NA))
colnames(cor_table_SST) = names(reg_df)
rownames(cor_table_SST) = ahba_probes$gene_symbol

# (pre-allocate) matrix to store correlation of every gene to PVALB
cor_table_PVALB           = data.frame(matrix(nrow=nrow(ahba_probes), ncol=length(reg_df), NA))
colnames(cor_table_PVALB) = names(reg_df)
rownames(cor_table_PVALB) = ahba_probes$gene_symbol

# results are tables with rows=20k probes, cols=8 regions
dim(cor_table_SST)
head(cor_table_PVALB)

# rotate z-transformed expression (mean expression of a gene across samples, m=0; sd=1)
ctx_micro_scale = as.data.frame(t(ctx_data_scale))


# give the user some info about the matrix
dim(ctx_micro_scale)
ctx_micro_scale[1:10,1:10]
mean(colMeans(ctx_micro_scale))

# "ahba_sorted" gives gene information (e.g. hgnc, entrez) for the expression matrix
hgnc_symbols = colnames(ctx_micro_scale)
ahba_sorted  = ahba_probes[match(hgnc_symbols, ahba_probes$gene_symbol),]
length(which(ahba_sorted$gene_symbol == hgnc_symbols))
head(ahba_sorted)


# cortical correlations of every gene to sst/pvalb
tmp_micro = ctx_micro_scale
sst_cor   = t(cor(tmp_micro[['SST']], tmp_micro))
pvalb_cor = t(cor(tmp_micro[['PVALB']], tmp_micro))


# relative position of SST/PVALB
length(which(sst_cor <= sst_cor[rownames(sst_cor) == 'PVALB',]))
length(which(pvalb_cor <= pvalb_cor[rownames(pvalb_cor) == 'SST',]))

# histogram of SST expression
fig_out = paste0(base_dir, '/figures/PaperFig_sst_expr_distribution.pdf')
fig_out
pdf(fig_out, width=3, height=2)
sst_expr_df = data.frame(sst_expr=sst_expr)
ggplot(sst_expr_df, aes(sst_expr)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('Cortical SST Expression')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=12),
          axis.text.y = element_text(colour="black", size=12),
          axis.title.x = element_text(colour="black", size=12),
          axis.title.y = element_text(colour="black", size=12)) +
    ylab("Density") +
    scale_x_continuous(expand = c(0, 0), limits=c(-5,5), breaks = seq(-5,5,2.5))
dev.off()


# histogram of PVALB expression
fig_out = paste0(base_dir, '/figures/PaperFig_pvalb_expr_distribution.pdf')
fig_out
pdf(fig_out, width=3, height=2)
pvalb_expr_df = data.frame(pvalb_expr=pvalb_expr)
ggplot(pvalb_expr_df, aes(pvalb_expr)) +
    geom_density(color="grey20", fill='lightblue', alpha=1) +
    ggtitle(paste0('Cortical PVALB Expression')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=12),
          axis.text.y = element_text(colour="black", size=12),
          axis.title.x = element_text(colour="black", size=12),
          axis.title.y = element_text(colour="black", size=12)) +
    ylab("Density") +
    scale_x_continuous(expand = c(0, 0), limits=c(-5,5), breaks = seq(-5,5,2.5))
dev.off()



# SST/PVALB corr for each of the 8 AHBA regions
pvarr = NULL
donor_probes = donorDat$probes

# we only want to analyze dorsal thalamus
thal_reg = c('DTA','ILc','LGd','DTLd','DTLv','MG','DTM','DTP','ILr')

for (reg in names(reg_df) ){
    write(reg,'')

    # mean-normalize expression seperately for each region/subject
    cur_reg_expr = NULL
    for (donor in donor_nums){
        write(donor,'')

        # region samples for this donor
        donor_reg_idxs = intersect(which(donorDat$samples$brain == donor), which(donorDat$samples$top_level == reg))

        # subset microarray and sample info
        donor_micro  = donorDat$micro[,donor_reg_idxs]
        donor_samps  = donorDat$samples[donor_reg_idxs,]

        # make sure we only analyze cortical samples that were within 4mm of the cortical sheet
        if (reg == 'CTX'){
            keep_samples = which(donor_samps$well_id %in% usable_ctx_samples$well_id)
            donor_samps  = donor_samps[keep_samples,]
            donor_micro  = donor_micro[,keep_samples]

        # only keep dorsal thalamic samples (defined above)
        } else if (reg == 'THAL'){
            keep_samples = which(donor_samps$structure_acronym %in% thal_reg)
            donor_samps  = donor_samps[keep_samples,]
            donor_micro  = donor_micro[,keep_samples]
        }

        reg_micro_scale = as.data.frame(apply(donor_micro, 1, scale, center=TRUE, scale=TRUE))
        colnames(reg_micro_scale) = donor_probes$gene_symbol
        cur_reg_expr = rbind(cur_reg_expr, reg_micro_scale)
    }
    cor(cur_reg_expr$SST, cur_reg_expr$PVALB)

    # correlation of sst/plvab
    print(dim(cur_reg_expr))
    corstats = cor.test(cur_reg_expr$SST, cur_reg_expr$PVALB)
    pvarr    = c(pvarr, corstats$p.value)
    
    # check stability with non-parametric stats
    print(cor.test(cur_reg_expr$SST, cur_reg_expr$PVALB))
    print(cor.test(cur_reg_expr$SST, cur_reg_expr$PVALB, method='spearman'))

    # Plot correlation of SST/PVALB for this region
    out_file = paste0(figure_dir, '/', reg, '_pvalb_vs_sst_expanded.pdf')
    plot_expression(cur_reg_expr, out_file, reg, reg_df)

    # Identify genes that are most correlated to SST / PVALB
    cor_table_SST[[reg]]   = t(cor(cur_reg_expr$SST, cur_reg_expr))
    cor_table_PVALB[[reg]] = t(cor(cur_reg_expr$PVALB, cur_reg_expr))
}
# reformat the region-wise data
sstpvalb_cors   = cor_table_SST[which(rownames(cor_table_SST) == 'PVALB'),]
sstpvalb_cors_t = data.frame(PVALB=t(sstpvalb_cors), region=colnames(sstpvalb_cors))
sstpvalb_cors_t = sstpvalb_cors_t[sstpvalb_cors_t$region != 'CTX',]
sstpvalb_cors_t = sstpvalb_cors_t[order(sstpvalb_cors_t$PVALB),]
sstpvalb_cors_t$region = factor(sstpvalb_cors_t$region, levels=sstpvalb_cors_t$region)
sstpvalb_cors_t
subcort_pval_df = data.frame(pval=pvarr[2:8], reg=names(reg_df))
subcort_pval_df$qval = p.adjust(subcort_pval_df$pval, method='BH')


# reorganize the dataframe for plotting
sstpvalb_cors_plot = sstpvalb_cors_t
colnames(sstpvalb_cors_plot) = c('cor','region')
sstpvalb_cors_plot = rbind(sstpvalb_cors_plot, sstpvalb_cors_plot)
sstpvalb_cors_plot$ref_gene = c(rep('PVALB',7), rep('SST',7))
sstpvalb_cors_plot$cor = as.numeric(as.character(sstpvalb_cors_plot$cor))


# region-wise distributions of all gene correlations to SST
cor_table_SST$gene = rownames(cor_table_SST)
cor_table_SST_long = gather(cor_table_SST, region, cor, CTX:GP)
cor_table_SST_long$ref_gene = 'SST'

# region-wise distributions of all gene correlations to PVALB
cor_table_PVALB$gene = rownames(cor_table_PVALB)
cor_table_PVALB_long = gather(cor_table_PVALB, region, cor, CTX:GP)
cor_table_PVALB_long$ref_gene = 'PVALB'

# combine PVALB/SST gene-wise correlations
cor_table_BOTH_long = rbind(cor_table_SST_long, cor_table_PVALB_long)
cor_table_BOTH_long = cor_table_BOTH_long %>% filter(region != 'CTX')
cor_table_BOTH_long = cor_table_BOTH_long %>% filter(cor != 1)

# convert region array to factor so they plot in the right order
sstpvalb_cors_t$region     = factor(sstpvalb_cors_t$region, levels=sstpvalb_cors_t$region)
cor_table_BOTH_long$region = factor(cor_table_BOTH_long$region, levels=sstpvalb_cors_t$region)

CairoPDF(paste0(figure_dir, '/SuppFig2_subcortical_Distributions_SSTvPVALB_correlations.pdf'), height=2.5, width=4)
ggplot(aes(y = cor, x = region, fill=ref_gene), data = cor_table_BOTH_long) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values=c('#5DA5DA', '#F15854')) +
            theme_minimal() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              axis.ticks.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x = element_blank(),
              axis.text.x = element_text(colour="black",size=16),
              axis.text.y = element_text(colour="black",size=10),
              axis.title.x = element_text(colour="black",size=10),
              axis.title.y = element_text(colour="black",size=10)) +
    expand_limits(y = c(-1, 1)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-1, 1, by = .5)) +
    geom_point(data=sstpvalb_cors_plot, aes(y = cor, x = region, fill=ref_gene), size=2, shape=21, stroke=1, color='black', fill='white')
dev.off()



# Collapse PVALB/SST relative expression by anatomically defined regions (e.g. Nucleus Accumbens, CA1, etc.)
plot_size_df = as.data.frame(matrix(NA, 2, 8))
colnames(plot_size_df) = names(reg_df)

# each region has different number of subregions, so we want to adjust figure dimensions for size uniformity
plot_size_df['CTX']  = rbind(c(6,3.5), c(4,-6))
plot_size_df['THAL'] = rbind(c(4,3), c(5,-5))
plot_size_df['HIPP'] = rbind(c(3,2), c(5,-5))
plot_size_df['STR']  = rbind(c(3,2.5), c(4.03,-4))
plot_size_df['VTA_SN'] = rbind(c(3,2.5), c(3,-3))
plot_size_df['AMYG'] = rbind(c(3,2.8), c(3,-3))
plot_size_df['HYPO'] = rbind(c(4,3), c(4,-4))
plot_size_df['GP']   = rbind(c(2.5,3), c(3,-3.05))


# Create sub-region PVALB/SST plots
donor_probes = donorDat$probes
all_region_data = NULL
# dorsal thalamus
thal_reg = c('DTA','ILc','LGd','DTLd','DTLv','MG','DTM','DTP','ILr')
for (reg in names(reg_df)[2:8]){
    print(reg)
    write(reg,'')

    # mean-normalize expression seperately for each region/subject
    cur_reg_expr = NULL
    cur_samples  = NULL
    for (donor in donor_nums){
        write(donor,'')
        donor_reg_idxs = intersect(which(donorDat$samples$brain == donor), which(donorDat$samples$top_level == reg))

        donor_micro  = donorDat$micro[,donor_reg_idxs]
        donor_samps  = donorDat$samples[donor_reg_idxs,]

        if (reg == 'CTX'){
            keep_samples = which(donor_samps$well_id %in% usable_ctx_samples$well_id)
            donor_samps  = donor_samps[keep_samples,]
            donor_micro  = donor_micro[,keep_samples]
        } else if (reg == 'THAL'){
            keep_samples = which(donor_samps$structure_acronym %in% thal_reg)
            donor_samps  = donor_samps[keep_samples,]
            donor_micro  = donor_micro[,keep_samples]
        }
        reg_micro_scale = as.data.frame(apply(donor_micro, 1, scale, center=TRUE, scale=TRUE))
        colnames(reg_micro_scale) = donor_probes$gene_symbol
        cur_reg_expr = rbind(cur_reg_expr, reg_micro_scale)
        cur_samples  = rbind(cur_samples, donor_samps)
    }

    cur_reg_expr$acro = cur_samples$structure_acronym

    # split substantia nigra samples into reticulata/compacta
    if (reg == 'VTA_SN'){
        new_names = as.character(cur_samples$structure_name)
        new_names[grep('reticulata', new_names)] = 'SN_pars_reticulata'
        new_names[grep('compacta', new_names)] = 'SN_pars_compacta'
        cur_samples$region_clean = new_names

        # t-test for SST/PVALB differences in the VTA/SNr
        # VTA
        print(t.test(cur_reg_expr$SST[cur_samples$structure_acronym == 'VTA'], cur_reg_expr$PVALB[cur_samples$structure_acronym == 'VTA'], paired=T))
        cur_reg_expr %>% filter(acro=='VTA') %>% summarize(PV_mean=mean(PVALB), PV_sd=sd(PVALB), SST_mean=mean(SST), SST_sd=sd(SST))

        # SNR
        print(t.test(cur_reg_expr$SST[cur_samples$structure_acronym == 'SNR'], cur_reg_expr$PVALB[cur_samples$structure_acronym == 'SNR'], paired=T))
        cur_reg_expr %>% filter(acro=='SNR') %>% summarize(PV_mean=mean(PVALB), PV_sd=sd(PVALB), SST_mean=mean(SST), SST_sd=sd(SST))
    }
    cur_reg_expr$region_clean = cur_samples$region_clean

    # only analyze regions with more than 3 samples
    examine_regs = names(table(cur_samples$region_clean))[which(table(cur_samples$region_clean) > 3)]
    gene_1 = 'SST'
    gene_2 = 'PVALB'

    # pre-allocate matrix to store relative difference information for each region
    subregion_diffs = as.data.frame(matrix(nrow=length(examine_regs), ncol=2, NA))
    print(paste0(reg, ': sst v pvalb: ', cor(cur_reg_expr[[gene_1]], cur_reg_expr[[gene_2]])))

    sst_summ_stats = cur_reg_expr %>%
                        filter(region_clean %in% examine_regs) %>%
                        group_by(region_clean) %>%
                        summarize(Median=median(SST),Min=min(SST),Max=max(SST),Q1=quantile(SST,.25),Q3=quantile(SST,.75),gene='SST',SSTmPVALB=median(SST)-median(PVALB))
    pvalb_summ_stats = cur_reg_expr %>%
                        filter(region_clean %in% examine_regs) %>%
                        group_by(region_clean) %>%
                        summarize(Median=median(PVALB),Min=min(PVALB),Max=max(PVALB),Q1=quantile(PVALB,.25),Q3=quantile(PVALB,.75),gene='PVALB',SSTmPVALB=median(SST)-median(PVALB))
    plot_data = rbind(sst_summ_stats, pvalb_summ_stats)
    all_region_data = rbind(all_region_data, plot_data)

    # sort regions by median sst/pvalb difference
    ordered_regs = as.character(unique(plot_data$region_clean[rev(order(plot_data$SSTmPVALB))]))
    plot_data$region_clean    = factor(plot_data$region_clean, levels=ordered_regs)

    width   = plot_size_df[[reg]][1]
    height  = plot_size_df[[reg]][1,2]
    plot_me = plot_data
    plot_me$gene = factor(plot_me$gene, levels=c('SST','PVALB'))
    
    CairoPDF(paste0(figure_dir, '/Paper_SuppFig_ahba_', reg, '_pvalb_vs_sst.pdf'), width, height)
    plot(ggplot(data=plot_me, aes(x=region_clean, y=Median, color=gene)) +
         geom_pointrange(aes(ymin=plot_me$Q1, ymax=plot_me$Q3), position=position_dodge(0.4), fatten = .75, size=.8) +
         geom_linerange(aes(ymin=plot_me$Min, ymax=plot_me$Max), position=position_dodge(0.4), size=.2) +
         theme_minimal() +
            scale_color_manual(values=c('#F15854','#5DA5DA', '#008000')) +
            theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.ticks.length=unit(.15, "cm"),
            axis.ticks = element_line(colour = "black"),
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size=20),
            axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
            axis.line.x = element_blank(),
            axis.text.x = element_text(colour="black",angle=90, hjust=1, size=5),
            axis.text.y = element_text(colour="black",size=10),
            axis.title.x = element_text(colour="black",size=10),
            axis.title.y = element_text(colour="black",size=10)) +
         scale_y_continuous(expand = c(0,0), limits=c(plot_size_df[[reg]][2,2], plot_size_df[[reg]][2,1])) +
         xlab("") +
         ylab("Relative SST-PVALB expression") +
         ggtitle(reg_df[[reg]]))
    dev.off()

    sst   = plot_me %>% filter(gene=='SST')
    pvalb = plot_me %>% filter(gene=='PVALB')
    cor.test(sst$Median, pvalb$Median)
}



# compare relative expression of SST and PVALB to 
# ground truth Rodent Cell Densities to AHBA SST-PVALB

# add AHBA acronym info to the "all_region_data" dataframe
acro_arr = NULL
for (reg in all_region_data$region_clean) {
    if (reg == 'SN_pars_compacta'){
        acro = 'SNc'
    } else if (reg == 'SN_pars_reticulata'){
        acro = 'SNr'
    } else {
        acro = donorDat$samples$structure_acronym[donorDat$samples$region_clean == reg][1]
    }
    acro_arr = c(acro_arr, acro)
}
all_region_data$acro = acro_arr



# read region-wise SST PVALB density estimates from Kim et al., 2018
kim_density    = read_csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/ref_files/kim_density.csv')
kim_density_df = kim_density[!is.na(kim_density$AHBA_acro),]

# average cell density by sub-area and collapse across sex
kim_density_df = kim_density_df %>% rowwise() %>%
                    mutate(Kim_PV_avg=mean(c(PV_M_mean, PV_F_mean), na.rm=T)) %>%
                    mutate(Kim_SST_avg=mean(c(SST_M_mean, SST_F_mean), na.rm=T))
kim_density_df = kim_density_df %>% group_by(AHBA_acro) %>% summarise(Kim_PV_avg=mean(Kim_PV_avg), Kim_SST_avg=mean(Kim_SST_avg))


# get AHBA region-wise Median expression in the same format
median_dat       = all_region_data[c('gene','Median','acro')] %>% spread(gene,Median)
caudopu_avg      = as.data.frame(t(colMeans(median_dat[median_dat$acro %in% c('BCd','HCd','TCd','Pu'),c('PVALB','SST')])))
caudopu_avg$acro = 'BCd;HCd;TCd;Pu'
ILM_avg          = as.data.frame(t(colMeans(median_dat[median_dat$acro %in% c('ILr','ILc'),c('PVALB','SST')])))
ILM_avg$acro     = 'Ilr;Ilc'
median_dat = rbind(median_dat, rbind(caudopu_avg, ILM_avg))


# combine rodent/AHBA data
merged = merge(x=kim_density_df, y=median_dat, by.x='AHBA_acro', by.y='acro')


# Difference and ratio of rodent SST and PVALB
merged$Kim_SST_PV   = merged$Kim_SST_avg - merged$Kim_PV_avg
merged$Kim_SST_PV_z = scale(merged$Kim_SST_avg) - scale(merged$Kim_PV_avg)
merged$AHBA_SST_PV  = merged$SST - merged$PVALB
merged$Kim_SST_PV_ratio = merged$Kim_SST_avg / merged$Kim_PV_avg


cor.test(merged$Kim_SST_PV, merged$AHBA_SST_PV, method='pearson')
cor.test(merged$Kim_SST_PV, merged$AHBA_SST_PV, method='spearman')
cor.test(merged$Kim_SST_PV_ratio, merged$AHBA_SST_PV, method='spearman')


# Figure 3b
CairoPDF('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_subcort_noscale.pdf', 4, 4)
p = ggplot(merged, aes(Kim_SST_PV, AHBA_SST_PV)) +
      geom_point(size=3, alpha=.5, stroke=0, colour='black') +
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
    xlab("Rodent SST-PVALB Cell Density (Kim et al.)") +
    ylab("AHBA SST(z) - PVALB(z) Expression") +
    geom_smooth(method='lm', se=FALSE, linetype = "dashed") +
    expand_limits(x = c(-10000, 30000), y = c(-5, 5)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(-10000, 30000, by = 10000)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-5, 5, by = 2.5))
  print(p)
dev.off()





