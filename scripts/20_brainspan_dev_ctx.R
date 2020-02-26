library(tidyverse)
library(WGCNA)
library(sva)
library(Cairo)

# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
data_dir = paste0(base_dir, '/data/brainspan')


# read expression data
gene_columns = read_csv(paste0(data_dir, '/columns_metadata.csv'))
gene_rows    = read_csv(paste0(data_dir, '/rows_metadata.csv'))
expr_matrix  = read_csv(paste0(data_dir, '/expression_matrix.csv'), col_names = FALSE)
expr_matrix  = expr_matrix %>% dplyr::select(-one_of('X1'))
gene_columns$id = 1:nrow(gene_columns)

# remove genes without an entrez id
expr_subset = expr_matrix[!is.na(gene_rows$entrez_id),]
rows_subset = gene_rows[!is.na(gene_rows$entrez_id),]


# log2 transform expression data and collapse duplicated transcripts
log2_expr        = log2(expr_subset + 1)
collapsed_matrix = collapseRows(log2_expr, rows_subset$entrez_id, rownames(log2_expr),
                                  method="MaxMean", connectivityBasedCollapsing=FALSE,
                                  connectivityPower=1, selectFewestMissing=TRUE, thresholdCombine=NA)
collapsed_expr   = collapsed_matrix$datETcollapsed
collapsed_expr   = as.data.frame(collapsed_expr)


# look-up dictionary to bin age ranges
age_bins = NULL
age_bins[['early_fetal']]    = c('8 pcw','9 pcw','12 pcw')
age_bins[['early_midfetal']] = c('13 pcw','16 pcw','17 pcw','19 pcw','21 pcw')
age_bins[['late_fetal']]     = c('24 pcw','25 pcw','26 pcw','35 pcw','37 pcw')
age_bins[['early_infancy']]  = c('4 mos')
age_bins[['late_infancy']]   = c('10 mos')
age_bins[['early_childhood']] = c('1 yrs','2 yrs','3 yrs','4 yrs')
age_bins[['mid_late_childhood']] = c('8 yrs', '11 yrs')
age_bins[['adolescence']]    = c('13 yrs', '15 yrs')
age_bins[['adult']]          = c('18 yrs','19 yrs','21 yrs','23 yrs','30 yrs','36 yrs','37 yrs','40 yrs')


# Cortical samples within each age bin
region_bins = NULL
region_bins[['early_fetal']]    = c("Ocx","M1C-S1C","STC","MFC","DFC","OFC","VFC","PCx","TCx","A1C","V1C","M1C","IPC","S1C")
region_bins[['early_midfetal']] = c("M1C","S1C","V1C","A1C","VFC","MFC","STC","OFC","IPC","DFC","ITC","M1C-S1C")
region_bins[['late_fetal']]     = c("DFC","S1C","A1C","IPC","OFC","STC","ITC","VFC","M1C","MFC","V1C")
region_bins[['early_infancy']]  = c("STC","V1C","ITC","MFC","M1C","DFC","OFC","A1C","VFC","S1C","IPC")
region_bins[['late_infancy']]   = c("S1C","IPC","STC","DFC","OFC","MFC","ITC","V1C")
region_bins[['early_childhood']] = c("M1C","DFC","STC","V1C","MFC","OFC","ITC","S1C","VFC","A1C","IPC")
region_bins[['mid_late_childhood']] = c("ITC","S1C","IPC","V1C","STC","VFC","M1C","DFC","A1C","MFC","OFC")
region_bins[['adolescence']] = c("MFC","IPC","VFC","A1C","ITC","STC","V1C","S1C","DFC","OFC","M1C")
region_bins[['adult']]       = c("M1C","V1C","A1C","ITC","VFC","MFC","STC","DFC","S1C","IPC","OFC")


# For each age bin:
# (1) combine data from different subjects,
# (2) remove donor effects if needed (combat)
# (3) linear regression, controlling for age
age_order    = NULL
dev_df       = NULL
sst_entrez   = '6750'
pvalb_entrez = '5816'
for ( bin in names(age_bins) ){
  write(bin, '')

  ages      = age_bins[[bin]]
  age_order = c(age_order, ages)
  ctx_acros = region_bins[[bin]]

  # tissue column info for cortical samples in this age group
  bin_col = gene_columns %>%
                rownames_to_column() %>%
                filter(age %in% ages) %>%
                filter(structure_acronym %in% ctx_acros)

  # expression info for this age bin
  bin_expr = collapsed_expr[bin_col$id]
   # get rid of genes that aren't expressed highly
  bin_expr = bin_expr[which(rowSums(bin_expr) > 2),]
  print(dim(bin_col))
  print(length(unique(bin_col$donor_id)))

  # only one subject in late-infancy, so don't run combat
  if (bin == 'late_infancy'){
    combat_expr = as.matrix(bin_expr)
  } else {
    combat_expr = ComBat(as.matrix(bin_expr), factor(bin_col$donor_name))
  }

  sst_expr   = scale(combat_expr[rownames(combat_expr) == sst_entrez,])
  pvalb_expr = scale(combat_expr[rownames(combat_expr) == pvalb_entrez,])

  bin_col$sst   = sst_expr
  bin_col$pvalb = pvalb_expr
  
  # these two bins only have one age value, so can't covary age
  if ( bin %in% c('early_infancy', 'late_infancy') ){
    out = summary(lm('sst_expr ~ pvalb_expr', bin_col))
  } else {
    out = summary(lm('sst_expr ~ pvalb_expr + age', bin_col))
  }
  coefs   = out$coefficients[rownames(out$coefficients) == 'pvalb_expr',]
  names(coefs) = c('est','se','tval','pval')
  cor_val = cor(sst_expr, pvalb_expr)
  dev_df  = rbind(dev_df, data.frame(t(coefs), bin))
}
dev_df

file_out = paste0(base_dir, '/figures/PaperFig_sst_pvalb_cor_brainspan.pdf')
CairoPDF(file_out)
ggplot() +
  geom_errorbar(data=dev_df, aes(x=bin, ymin=est-se, ymax=est+se), width=0.2) +
  geom_point(data=dev_df, aes(x=bin, y=est), size=4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    scale_y_continuous(limits=c(-1,.5), expand=c(0,0))
dev.off()







