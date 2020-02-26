library(tidyverse)
library(dplyr)
library(Cairo)

# function to map ontology IDs to larger overall categories in our analysis (e.g. CTX, HIPP)
find_top_level = function(ref_df, struct_id, ontology) {
  ontology_row  = which(ontology$id == struct_id )
  if (length(ontology_row) > 1){ # print feedback if more than one match, this should never be the case
    write("ERROR")
  }
  ontology_info = ontology[ontology_row,]
  splits    = strsplit(as.character(ontology_info$structure_id_path), '/')[[1]]
  reg.match = ref_df$id[which(top.level$id %in% splits)]
  if (length(reg.match) == 1){
    region_cat  = as.character(top.level$top_level[which(top.level$id %in% splits)])
    region_name = top.level[top.level$id == reg.match,]$name
    return(c(reg.match, region_cat, region_name))
  } else {
    return(c(NA,NA,NA))
  }
}


# project directories
interneuron_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
rnaseq_dir      = paste0(interneuron_dir, '/data/ahba/rna_seq')


# list with subject info
donors     = c('9861', '10021')
rna_seq_df = NULL
for (donor in donors){
    print(donor)
    rna_seq_df[[donor]] = NULL

    # donor specific directory
    donor_dir  = paste0(rnaseq_dir, '/rnaseq_donor',donor)

    # read sample annotations
    annot_dat  = read_csv(paste0(donor_dir, '/SampleAnnot.csv'))
    rna_seq_df[[donor]][['annot']] = annot_dat

    # subset to cortical only sample annotations
    ctx_idxs   = rna_seq_df[[donor]][['annot']][['main_structure']] %in% c('FL','OL','TL','PL')
    rna_seq_df[[donor]][['samples_ctx']] = rna_seq_df[[donor]][['annot']][ctx_idxs,]

    #
    ont_dat = read_csv(paste0(donor_dir, '/Ontology.csv'))
    rna_seq_df[[donor]][['ontology']] = ont_dat

    #
    gene_dat = read_csv(paste0(donor_dir, '/Genes.csv'))
    rna_seq_df[[donor]][['genes']] = gene_dat

    # expression data (transcripts per million)
    tpm_dat = as.data.frame(read_csv(paste0(donor_dir, '/RNAseqTPM.csv'), col_names=F))
    rownames(tpm_dat) = tpm_dat$X1
    tpm_dat$X1 = NULL

    rna_seq_df[[donor]][['tpm']] = tpm_dat
    rna_seq_df[[donor]][['tpm_ctx']] = tpm_dat[ctx_idxs]
    ctx_expr = rna_seq_df[[donor]][['tpm_ctx']]

    rna_seq_df[[donor]][['tpm_ctx_scale']] = as.data.frame(t(apply(ctx_expr, 1, scale, center=TRUE, scale=TRUE)))
}

# combine the two subjects experssion/sample data
all_ctx_expr = cbind(rna_seq_df[['9861']][['tpm_ctx']], rna_seq_df[['10021']][['tpm_ctx']])
all_ctx_scaleexpr = cbind(rna_seq_df[['9861']][['tpm_ctx_scale']], rna_seq_df[['10021']][['tpm_ctx_scale']])
all_ctx_samp = rbind(rna_seq_df[['9861']][['samples_ctx']], rna_seq_df[['10021']][['samples_ctx']])
all_ctx_samp = as.data.frame(all_ctx_samp)


# 184 cortex samples
n_samps = ncol(all_ctx_scaleexpr)

# for each gene, percentage of cortical samples with zero expression
zero_expr_cts  = rowSums(all_ctx_expr == 0)
zero_expr_perc = rowSums(all_ctx_expr == 0)/n_samps

# genes that are not expressed in at least 50% of cortical samples
low_expr_genes = which(zero_expr_perc <= .5)

# remove the low expression geenes
genes_filter = gene_dat[low_expr_genes,]
all_ctx_scaleexpr_filter = all_ctx_scaleexpr[low_expr_genes,]


# calc correlation between SST/PVALB
sst_idx    = which(genes_filter$gene_symbol == 'SST')
pvalb_idx  = which(genes_filter$gene_symbol == 'PVALB')
sst_expr   = as.numeric(all_ctx_scaleexpr_filter[sst_idx,])
pvalb_expr = as.numeric(all_ctx_scaleexpr_filter[pvalb_idx,])

cor.test(sst_expr, pvalb_expr, method='spearman')
sst_pvalb_cor_sp = cor.test(sst_expr, pvalb_expr, method='spearman')

cor.test(sst_expr, pvalb_expr, method='pearson')
sst_pvalb_cor = cor.test(sst_expr, pvalb_expr, method='pearson')

# organize df for plotting
plot_me  = data.frame(sst=sst_expr, pvalb=pvalb_expr, reg = all_ctx_samp$ontology_structure_acronym)

out_file = paste0(interneuron_dir, '/figures/PaperFig_CorrPlot_TPM_pvalb_sst_cor.pdf')
CairoPDF(out_file, width=4, height=4, fonts='Arial')
p = ggplot(plot_me, aes(sst, pvalb)) +
            geom_point(size=2, alpha=1, stroke=.5, colour='black', shape=21) +
            ggtitle('RNASeq Cortical Expression') +
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
            ylab("Normalized PVALB Expression (z)") +
            xlab("Normalized SST Expression (z)") +
            geom_smooth(method='lm', linetype = "dashed", fullrange = T, color='black') +
            expand_limits(x = c(-4,3), y = c(-2,4)) +
            scale_x_continuous(expand = c(0, 0), breaks = seq(-4,3, by = 1)) +
            scale_y_continuous(expand = c(0, 0), breaks = seq(-2,4, by = 1)) +
            coord_fixed()
print(p)
dev.off()
out_file


# summarize expression by each sub-region
regions_w_dat = rna_seq_df[[donor]][['ontology']][rna_seq_df[[donor]][['ontology']]$id %in% rna_seq_df[[donor]][['annot']]$ontology_structure_id,]

sst_dat = plot_me %>%
            group_by(reg) %>%
            summarise(mean = mean(sst),
                      median = median(sst),
                      Q1 = quantile(sst,.25),
                      Q3 = quantile(sst,.75),
                      min = min(sst),
                      max = max(sst))
sst_dat$gene = 'sst'

# plot with the full name for each area, rather than the acryonym
full_names = NULL
for (r in sst_dat$reg){
    reg_names = ont_dat$name[ont_dat$acronym == r]
    reg_names = gsub('left|right|[,]','',reg_names)
    reg_names = gsub(' ','_',reg_names)
    full_names = c(full_names, reg_names[1])
    print('')
    print(unique(reg_names))
}
sst_dat$full_names = full_names
sst_dat$full_names = gsub('_',' ', gsub('__','_',as.character(sst_dat$full_names)))


pvalb_dat = plot_me %>%
            group_by(reg) %>%
            summarise(mean = mean(pvalb),
                      median = median(pvalb),
                      Q1 = quantile(pvalb,.25),
                      Q3 = quantile(pvalb,.75),
                      min = min(pvalb),
                      max = max(pvalb))
pvalb_dat$gene = 'pvalb'


# make an array of full region names rather than acronyms
full_names = NULL
for (r in sst_dat$reg){
    reg_names = ont_dat$name[ont_dat$acronym == r]
    reg_names = gsub('left|right|[,]','',reg_names)
    reg_names = gsub(' ','_',reg_names)
    full_names = c(full_names, reg_names[1])
    print('')
    print(unique(reg_names))
}
pvalb_dat$full_names = full_names
pvalb_dat$full_names = gsub('_',' ', gsub('__','_',as.character(pvalb_dat$full_names)))


# region-wise difference in SST and PVALB medians
diffs          = sst_dat$median - pvalb_dat$median
names(diffs)   = sst_dat$full_names
plot_reg_order = names(diffs)[rev(order(diffs))]


# get ready for plotting
plot_me = rbind(sst_dat, pvalb_dat)
plot_me$full_names = factor(as.character(plot_me$full_names), levels = plot_reg_order)
sst_dat$reg == pvalb_dat$reg
cor.test(sst_dat$median, pvalb_dat$median)
cor.test(sst_dat$median, pvalb_dat$median, method='spearman')

# clean up summary stat data frame
plot_me$gene = factor(plot_me$gene, levels=c('sst','pvalb'))
out_file = paste0(interneuron_dir, '/figures/PaperFig_RegionWise_TPM_pvalb_sst_cor.pdf')
CairoPDF(out_file, width=6, height=5)
p = ggplot(data=plot_me, aes(x=full_names, y=median, color=gene)) +
     geom_pointrange(aes(ymin=plot_me$Q1, ymax=plot_me$Q3), position=position_dodge(0.4), size=.8) +
     geom_linerange(aes(ymin=plot_me$min, ymax=plot_me$max), position=position_dodge(0.4), size=.2) +
     theme_minimal() + scale_color_manual(values=c('#F15854','#5DA5DA', '#008000')) +
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
     scale_y_continuous(expand = c(0,0), limits=c(-4,4)) +
     xlab("") +
     ylab("Relative SST-PVALB expression")
print(p)
dev.off()
out_file



# benchmark SST and PVALB correlation against all other genes
sst_expr   = as.numeric(all_ctx_scaleexpr_filter[sst_idx,])
pvalb_expr = as.numeric(all_ctx_scaleexpr_filter[pvalb_idx,])


# all genes to cortical SST
sst_cors      = as.data.frame(cor(sst_expr, t(all_ctx_scaleexpr_filter), method='pearson'))
sst_cors_sort = sst_cors[rev(order(sst_cors))]
head(sst_cors_sort)

# get rid of SSTs correlatino to itself
sst_cors_sort = sst_cors_sort[names(sst_cors_sort) != 'SST']
plot_df       = as.data.frame(sst_cors_sort)

# proportion of genes that are more negative than pvalb
length(which(plot_df <= as.numeric(plot_df[names(plot_df) == 'PVALB'])) / length(plot_df))

# get histogram ready
plot_me = data.frame(genes=names(plot_df), cor=as.numeric(plot_df))
sst_mean = mean(plot_me$cor)


# plot hist!
out_file = paste0(interneuron_dir, '/figures/PaperFig_TPM_SSTcor_distr.pdf')
CairoPDF(out_file, width=3, height=2)
p = ggplot(plot_me, aes(cor)) +
    geom_density(color="grey20", fill='#F15854', alpha=1) +
    ggtitle('sst') +
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
    xlab("Avg Cross-Region Correlation to SST") +
    expand_limits(y=c(0,2)) +
    scale_x_continuous(expand = c(0, 0), limits=c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = sst_pvalb_cor$estimate, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = sst_mean, color='black', size=.5, linetype='longdash') +
      annotate("text", x=sst_pvalb_cor$estimate, y=1, label='PVALB', color ='#5DA5DA')
print(p)
dev.off()
out_file





# all genes correlated to PVALB
pvalb_cors = as.data.frame(cor(pvalb_expr, t(all_ctx_scaleexpr_filter), method='pearson'))
pvalb_cors_sort = pvalb_cors[rev(order(pvalb_cors))]

# get rid of PVALbs correlation to itself
pvalb_cors_sort = pvalb_cors_sort[names(pvalb_cors_sort) != 'PVALB']
plot_df = as.data.frame(pvalb_cors_sort)

# proportion of genes that are more negative than pvalb
length(which(plot_df <= as.numeric(plot_df[names(plot_df) == 'SST']))) / length(plot_df)

# get ready for hist
plot_me    = data.frame(genes=names(plot_df), cor=as.numeric(plot_df))
pvalb_mean = mean(plot_me$cor)

out_file = paste0(interneuron_dir, '/figures/PaperFig_TPM_PVALBcor_distr.pdf')
CairoPDF(out_file, width=3, height=2)
p = ggplot(plot_me, aes(cor)) +
    geom_density(color="grey20", fill='#5DA5DA', alpha=1) +
    ggtitle('sst') +
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
    xlab("Avg Cross-Region Correlation to PVALB") +
    expand_limits(y=c(0,2)) +
    scale_x_continuous(expand = c(0, 0), limits=c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by = 1)) +
    geom_vline(xintercept = sst_pvalb_cor$estimate, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pvalb_mean, color='black', size=.5, linetype='longdash') +
      annotate("text", x=sst_pvalb_cor$estimate, y=1, label='SST', color ='#F15854')
print(p)
dev.off()
out_file



# compare SST and PVALB correlation to all other parwise correlations
# all pairwise correlations
all_cors = cor(t(all_ctx_scaleexpr_filter))
x = all_cors[upper.tri(all_cors)]
all_cors = data.frame(cor=x)

overall_mean = mean(all_cors$cor)
mean(all_cors$cor)
sd(all_cors$cor)

length(which(all_cors$cor <= sst_pvalb_cor$estimate))/length(all_cors$cor)

out_file = paste0(interneuron_dir, '/figures/PaperFig_TPM_AllpairwiseCorrs.pdf')
CairoPDF(out_file, width=3, height=2)
p = ggplot(all_cors, aes(cor)) +
    geom_density(color="grey20", fill='#5DA5DA', alpha=1) +
    ggtitle('sst') +
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
    xlab("All pairwise correlations") +
    expand_limits(y=c(0,2)) +
    scale_x_continuous(expand = c(0, 0), limits=c(-1, 1), breaks = c(-1, -.5, 0, .5, 1)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2, by=1)) +
    geom_vline(xintercept = sst_pvalb_cor$estimate, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = overall_mean, color='black', size=.5, linetype='longdash') +
      annotate("text", x=sst_pvalb_cor$estimate, y=1, label='SST_v_PVALB', color ='#F15854')
print(p)
dev.off()
out_file






