library(ggplot2)
library(affy)
library(WGCNA)
library(GEOquery)
library(sva)
library(tidyverse)
library(Cairo)


# Bernard Cortical Data
base_dir    = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
bernard_dir = paste(base_dir, '/data/bernard', sep = '')


# RMA and batch normalized log10 expression data is available on NIH GSE Data website
bernard_gse = getGEO(filename=paste0(bernard_dir, '/GSE31613_series_matrix.txt'), GSEMatrix=TRUE)


# extract the expression values, convert to log2
bernard_expr_matrix = bernard_gse@assayData$exprs
dim(bernard_expr_matrix)
max(bernard_expr_matrix)
bernard_expr_log2 = log2(10^bernard_expr_matrix)
max(bernard_expr_log2)
colnames(bernard_expr_log2) = bernard_gse@phenoData@data$source_name_ch1 # set region names to column names


# subset to the 4 primates with complete sampling
#   99C039 - only has middle temporal area sampling
#   01D150 - only has OFC, ACC, and DLPFC sampling
exp_animals     = c('018G', '98C096', '99C018', 'CR5C1')
exp_animal_idxs = sort(unlist(lapply(exp_animals, grep, x=as.character(bernard_gse@phenoData@data$title))))

bernard_expr_use = bernard_expr_log2[,exp_animal_idxs]
bernard_expr_use[1:10,1:10]
pheno_data_matrix = bernard_gse@phenoData@data[exp_animal_idxs,]
head(pheno_data_matrix)


# organize tissue sample information
pheno_df = NULL
pheno_df$title       = as.character(pheno_data_matrix$title)
pheno_df$tissue_type = gsub('tissue: ', '', as.character(pheno_data_matrix$characteristics_ch1.1))
pheno_df$layer       = gsub('tissue region: ', '', as.character(pheno_data_matrix$characteristics_ch1.2))
pheno_df$gender      = pheno_data_matrix[['gender:ch1']]
pheno_df = as_tibble(pheno_df)

# do a bit of renaming
cort_regions = c('Anterior cingulate gyrus','Orbitofrontal cortex', 'Primary auditory ctx',
                      'Dorsolateral Prefrontal Cortex','Primary somatosensory cortex',
                      'Primary motor cortex', 'Temporal area', 'Middle temporal area','V1', 'V2')
pheno_df$tissue_type = gsub('middle', 'Middle', pheno_df$tissue_type)


# subset data to only include cortical regions
cort_idxs = sort(unlist(lapply(cort_regions, grep, x=pheno_df$tissue_type)))
pheno_df  = pheno_df[cort_idxs,]
cort_expr = bernard_expr_use[,cort_idxs]

# make donor array
pheno_df$donor = NA
for (d in exp_animals){
  pheno_df$donor[grep(d, pheno_df$title)] <- d
}


# Only analyze high probability human/macaque gene homologs identified by Bakken et al., 2016
bakken_df = read.csv(paste(base_dir, '/data/bakken/bakken_probe_mapping.csv', sep = ''))
bakken_df = bakken_df %>% filter(!is.na(human_entrezid))# bakken probes with human entrez
bakken_df = bakken_df %>% filter(keep_for_analysis==TRUE)


# add human gene names to the Bakken probe data frame
bakken_order    = bakken_df[order(bakken_df$probeid),]
cort_subset     = cort_expr[rownames(cort_expr) %in% bakken_df$probeid,]
cort_expr_order = cort_subset[order(rownames(cort_subset)),]
length(rownames(cort_expr_order) == bakken_order$probeid)

# average probes if they have non-unique mapping to a human gene (only affects a handful of probes)
expr_collapse_tmp  = collapseRows(cort_expr_order, as.character(bakken_order$human_genesymbol), rownames(cort_expr_order))
expr_collapse      = expr_collapse_tmp$datETcollapsed
expr_combat        = ComBat(expr_collapse, pheno_df$donor, prior.plots=F)
expr_combat_scale  = as.data.frame(apply(expr_combat, 1, scale, center=TRUE, scale=TRUE))

# average expression by region/donor
avg_mat = NULL
for (donor in unique(pheno_df$donor)){
  for (reg in cort_regions){
    don_idxs = which(pheno_df$donor == donor)
    reg_idxs = grep(reg, pheno_df$tissue_type)
    match_idxs = intersect(don_idxs, reg_idxs)
    if ( length(match_idxs) > 0 ){
      name = paste(donor, reg, sep='_')
      write(name,'')
      avg_mat[[name]] = colMeans(expr_combat_scale[match_idxs, ])
    }
  }
}
avg_mat = as.data.frame(avg_mat)

# cortical PVALB to SST correlation
cor.test(as.numeric(avg_mat[rownames(avg_mat) == 'PVALB',]), as.numeric(avg_mat[rownames(avg_mat) == 'SST',]))
cor.test(as.numeric(avg_mat[rownames(avg_mat) == 'PVALB',]), as.numeric(avg_mat[rownames(avg_mat) == 'SST',]), method='spearman')

pvalb.expr = avg_mat[which(rownames(avg_mat)== 'PVALB'),]
sst.expr   = avg_mat[which(rownames(avg_mat)== 'SST'),]
cor(t(pvalb.expr), t(sst.expr))

# Plot correlation
avg_mat_t  = as.data.frame(t(avg_mat))
figure_dir = paste0(base_dir, '/figures')
CairoPDF(paste0(figure_dir, '/PaperFig_primate_CTX_pvalb_sst_overall_corr.pdf'))
ggplot(avg_mat_t, aes(PVALB, SST)) + geom_point(size=5, stroke=0, color='#2EA121') +
  ggtitle('CTX') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
        axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
        axis.text.x = element_text(colour="grey20",size=20),
        axis.text.y = element_text(colour="grey20",size=20),
        axis.title.x = element_text(colour="grey20",size=20),
        axis.title.y = element_text(colour="grey20",size=20)) +
  ylab("Normalized PVALB Expression (z)") +
  xlab("Normalized SST Expression (z)") +
  expand_limits(x = c(-2, 2), y = c(-2, 2)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-2, 0, 2)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(-2, 0, 2)) +
  geom_smooth(method='lm', se=FALSE)
dev.off()
cor.test(avg_mat_t$PVALB, avg_mat_t$SST)


# summary stats of SST/PVALB for each region
region_expr = NULL
for (reg in gsub(' ', '.', cort_regions)){
  reg_expr   = avg_mat[grep(reg, colnames(avg_mat))]
  sst_summ   = summary(as.numeric(t(reg_expr[rownames(reg_expr) == 'SST',])))
  pvalb_summ = summary(as.numeric(t(reg_expr[rownames(reg_expr) == 'PVALB',])))
  reg_df     = as.data.frame(rbind(sst_summ, pvalb_summ))
  reg_df$gene = c('SST','PVALB')
  reg_df$region = gsub('[.]',' ',reg)
  region_expr   = rbind(region_expr, reg_df)
}
rownames(region_expr) = NULL
colnames(region_expr) = c('Min','Q1','Median','Mean','Q3','Max','gene','region')


# plot!
CairoPDF(paste0(figure_dir, '/primate_CTX_pvalb_sst_circlboxplot_corr.pdf'))
sst_expr           = region_expr %>% filter(gene=='SST')
pvalb_expr         = region_expr %>% filter(gene=='PVALB')
ordered_regions    = sst_expr$region[rev(order(sst_expr$Median - pvalb_expr$Median))]
region_expr$region = factor(region_expr$region, levels=ordered_regions)
ggplot(data=region_expr, aes(x=region, y=Median, color=gene)) +
        geom_pointrange(aes(ymin=region_expr$Q1, ymax=region_expr$Q3), position=position_dodge(0.4), size=.8) +
        geom_linerange(aes(ymin=region_expr$Min, ymax=region_expr$Max), position=position_dodge(0.4), size=.2) +
        theme_minimal() + scale_color_manual(values=c('#5DA5DA','#F15854')) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              axis.ticks.x = element_blank(),
              plot.title = element_text(hjust = 0.5, size=20),
              #axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x = element_blank(),
              axis.text.x = element_text(colour="black",angle=45, hjust=1, size=5),
              axis.text.y = element_text(colour="black",size=10),
              axis.title.x = element_text(colour="black",size=10),
              axis.title.y = element_text(colour="black",size=10)) +
        scale_y_continuous(expand = c(0,0), limits=c(-2,2)) +
        xlab("") +
        ylab("Relative SST-PVALB expression") +
        ggtitle('Cortex')
dev.off()



