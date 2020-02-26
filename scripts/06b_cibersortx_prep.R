library(SingleCellExperiment)
library(Seurat)
library(mclust)
library(dplyr)
library(GEOquery)
library(EWCE)
library(tidyverse)
library(biomaRt)


# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
umi_dir  = file.path(base_dir, 'data/singlecell/')

region  = 'VisualCortex'
in_file = paste0(umi_dir, 'lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed.Rdata')
in_file
load(verbose=T, file = in_file)
#Loading objects:
#  lake_vis_seuset
#  lake_vis_gene_dict

# load AHBA data
load(file=paste0(base_dir, '/data/ahba/ahba_data_object.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)

donor_nums = names(ahba_data)

# grab the probe information
ahba_probes = donorDat$probes

# ENTREZ IDs that match between the two datasets
id_keep = intersect(lake_vis_gene_dict$entrez_id, ahba_probes$entrez_id)

# rename the ahba gene name array
lake_vis_ahba_combo = merge(x=ahba_probes, lake_vis_gene_dict, by.x='entrez_id', by.y='entrez_id')


# prep SINGLE CELL SIGNATURE FILE
lake_vis_keep_idxs = which(lake_vis_gene_dict$entrez_id %in% id_keep)
lake_vis_micro = 10^as.matrix(lake_vis_seuset@data[lake_vis_keep_idxs,])
lake_vis_micro_orig = as.matrix(lake_vis_seuset@data[lake_vis_keep_idxs,])
lake_vis_scale = as.matrix(lake_vis_seuset@scale.data[lake_vis_keep_idxs,])

# replace the lake gene names with the ones that have been aligned to the AHBA data
lake_vis_subset_df = data.frame(genes=rownames(lake_vis_micro), lakeidx=1:nrow(lake_vis_micro))
lake_vis_subset_merge_df = merge(x=lake_vis_subset_df, y=lake_vis_ahba_combo, by.x='genes', by.y='lake_orig_Symbol')
lake_vis_subset_merge_df = lake_vis_subset_merge_df[order(lake_vis_subset_merge_df$lakeidx),]

# get rid of some duplicates
keep_idxs2     = which(!duplicated(lake_vis_subset_merge_df$gene_symbol))
lake_vis_micro = lake_vis_micro[keep_idxs2,]
lake_vis_micro_orig = lake_vis_micro_orig[keep_idxs2,]
lake_vis_scale = lake_vis_scale[keep_idxs2,]
lake_vis_probe_info = lake_vis_subset_merge_df[keep_idxs2,]

# write data
save(lake_vis_micro_orig, lake_vis_scale, lake_vis_probe_info, file=paste0(umi_dir, 'lake2018/lake_vis_ahba_match.Rdata'))
load(file=paste0(umi_dir, 'lake2018/lake_vis_ahba_match.Rdata'), verbose=T)



# CIBERSORT formatting
lake_vis_micro = 10^lake_vis_micro_orig

# CIBERSORT formatting
vis_dexmat = cbind(lake_vis_probe_info$gene_symbol, lake_vis_micro)
cell_ident = unlist(lapply(colnames(lake_vis_seuset@data), function(x) strsplit(x,'_')[[1]][[1]]))

colnames(vis_dexmat) = c('!Sample_title', cell_ident)
sample_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_orig_celllabels.txt')
write.table(x=vis_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')

cell_ident = as.character(lake_vis_seuset@ident)
cell_ident[grep('Ex6',cell_ident)] = 'Ex6'
cell_ident[grep('Ex3',cell_ident)] = 'Ex3'
cell_ident[grep('Ex5',cell_ident)] = 'Ex5'
cell_ident[grep('In1',cell_ident)] = 'In1'
cell_ident[grep('In4',cell_ident)] = 'In4'
cell_ident[grep('In6',cell_ident)] = 'PVALB'
cell_ident[grep('In7|In8',cell_ident)] = 'SST'

colnames(vis_dexmat) = c('!Sample_title', cell_ident)

# write
sample_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample.txt')
write.table(x=vis_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')


# CIBERSORT phenotyping data
row = 1
pheno_mat = matrix(2, nrow=length(unique(unique(lake_vis_seuset@ident))), ncol=length(lake_vis_seuset@ident))
for (cell in unique(lake_vis_seuset@ident)){
    pheno_mat[row,lake_vis_seuset@ident == cell] = 1
    row = row + 1
}
uniq_cell_arr = unique(lake_vis_seuset@ident)
out_pheno_mat = cbind(as.character(unique(lake_vis_seuset@ident)), pheno_mat)
pheno_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_pheno.txt')
write.table(x=out_pheno_mat, file=pheno_file, quote=FALSE, row.names = FALSE, sep='\t', col.names=F)


# subset AHBA data and WRITE MIXTURE FILE
for (donor in donor_nums){
    write(donor, '')
    keep_idxs = which(donorDat$probes$gene_symbol %in% lake_vis_probe_info$gene_symbol)

    donor_ctx_idxs = intersect(which(donorDat$samples$brain == donor), which(donorDat$samples$top_level == 'CTX'))

    ahba_ctx_dat    = 2^donorDat$micro[keep_idxs,donor_ctx_idxs]
    ahba_ctx_probes = donorDat$probes[keep_idxs,]
    ahba_ctx_samples = donorDat$samples[donor_ctx_idxs,]

    ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
    colnames(ahba_ctx_out) = c('!Sample_title', ahba_ctx_samples$well_id)

    ciber_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/cibersortX'
    mixture_file = paste0(ciber_dir, '/', region,'_ahba_',donor,'_mixture.txt')
    write(mixture_file,'')
    write.table(x=ahba_ctx_out, file=mixture_file, quote=FALSE, row.names = FALSE, sep='\t')
}





# set up directories
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
umi_dir  = file.path(base_dir, 'data/singlecell/')

region = 'FrontalCortex'
load(verbose=T, file = paste0(umi_dir, 'lake2018/GSE97930_',region,'_snDrop-seq_UMI_Count_Matrix_08-01-2017_seurat_processed.Rdata'))
# Loading objects:
#  lake_dfc_seuset
#  lake_dfc_gene_dict
# lake_dfc_gene_dict =lake_gene_dict

# grab the probe information
ahba_probes = donorDat$probes

# ENTREZ IDs that match between the two datasets
id_keep = intersect(lake_dfc_gene_dict$entrez_id, ahba_probes$entrez_id)

# rename the ahba gene name array
#ahba_probes$ahba_hg37_GeneNames = ahba_probes$hg37_GeneNames
lake_dfc_ahba_combo = merge(x=ahba_probes, lake_dfc_gene_dict, by.x='entrez_id', by.y='entrez_id')


# prep SINGLE CELL SIGNATURE FILE
lake_dfc_keep_idxs = which(lake_dfc_gene_dict$entrez_id %in% id_keep)
lake_dfc_micro = 10^as.matrix(lake_dfc_seuset@data[lake_dfc_keep_idxs,])
lake_dfc_micro_orig = as.matrix(lake_dfc_seuset@data[lake_dfc_keep_idxs,])
lake_dfc_scale = as.matrix(lake_dfc_seuset@scale.data[lake_dfc_keep_idxs,])


# replace the lake gene names with the ones that have been aligned to the AHBA data
lake_dfc_subset_df = data.frame(genes=rownames(lake_dfc_micro), lakeidx=1:nrow(lake_dfc_micro))
lake_dfc_subset_merge_df = merge(x=lake_dfc_subset_df, y=lake_dfc_ahba_combo, by.x='genes', by.y='lake_orig_Symbol')
lake_dfc_subset_merge_df = lake_dfc_subset_merge_df[order(lake_dfc_subset_merge_df$lakeidx),]

# get rid of some duplicates
keep_idxs2     = which(!duplicated(lake_dfc_subset_merge_df$gene_symbol))
lake_dfc_micro = lake_dfc_micro[keep_idxs2,]
lake_dfc_micro_orig = lake_dfc_micro_orig[keep_idxs2,]
lake_dfc_scale = lake_dfc_scale[keep_idxs2,]
lake_dfc_probe_info = lake_dfc_subset_merge_df[keep_idxs2,]


# replace the lake gene names with the ones that have been aligned to the AHBA data
lake_dfc_subset_df = data.frame(genes=rownames(lake_dfc_micro), lakeidx=1:nrow(lake_dfc_micro))
lake_dfc_subset_merge_df = merge(x=lake_dfc_subset_df, y=lake_dfc_ahba_combo, by.x='genes', by.y='lake_orig_Symbol')
lake_dfc_subset_merge_df = lake_dfc_subset_merge_df[order(lake_dfc_subset_merge_df$lakeidx),]

# write data
save(lake_dfc_micro_orig, lake_dfc_scale, lake_dfc_subset_merge_df, file=paste0(umi_dir, 'lake2018/lake_dfc_ahba_match.Rdata'))
load(file=paste0(umi_dir, 'lake2018/lake_dfc_ahba_match.Rdata'), verbose=T)



# CIBERSORT formatting
#lake_vis_micro = 10^lake_vis_micro_orig
#lake_vis_micro = round(lake_vis_micro,6)
dfc_dexmat = cbind(lake_dfc_subset_merge_df$gene_symbol, lake_dfc_micro)
#cell_ident = as.character(lake_dfc_seuset@ident)
cell_ident = unlist(lapply(colnames(lake_dfc_seuset@data), function(x) strsplit(x,'_')[[1]][[1]]))

colnames(dfc_dexmat) = c('!Sample_title', cell_ident)
sample_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_orig_celllabels.txt')
write.table(x=dfc_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')

cell_ident[grep('In6', cell_ident)] = 'PVALB'
cell_ident[grep('In7|In8', cell_ident)] = 'SST'

colnames(dfc_dexmat) = c('!Sample_title', cell_ident)
sample_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample_orig_sstpvalb_combine.txt')
write.table(x=dfc_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')


cell_ident[grep('Ex6',cell_ident)] = 'Ex6'
cell_ident[grep('Ex3',cell_ident)] = 'Ex3'
cell_ident[grep('Ex5',cell_ident)] = 'Ex5'
cell_ident[grep('In1',cell_ident)] = 'In1'
cell_ident[grep('In4',cell_ident)] = 'In4'
#cell_ident[grep('In6',cell_ident)] = 'PVALB'
#cell_ident[grep('In7|In8',cell_ident)] = 'SST'

colnames(dfc_dexmat) = c('!Sample_title', cell_ident)

# write
sample_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_sample.txt')
write.table(x=dfc_dexmat, file=sample_file, quote=FALSE, row.names = FALSE, sep='\t')

# CIBERSORT phenotyping data
row = 1
pheno_matdf = matrix(2, nrow=length(unique(unique(lake_vis_seuset@ident))), ncol=length(lake_vis_seuset@ident))
for (cell in unique(lake_vis_seuset@ident)){
    pheno_mat[row,lake_vis_seuset@ident == cell] = 1
    row = row + 1
}
uniq_cell_arr = unique(lake_vis_seuset@ident)
out_pheno_mat = cbind(as.character(unique(lake_vis_seuset@ident)), pheno_mat)
pheno_file = paste0(base_dir, '/data/cibersortX/Lake_',region,'_ahba_matched_sc_pheno.txt')
write.table(x=out_pheno_mat, file=pheno_file, quote=FALSE, row.names = FALSE, sep='\t', col.names=F)



# subset AHBA data and WRITE MIXTURE FILE
for (donor in donor_nums){
    write(donor, '')
    keep_idxs = which(donorDat$probes$gene_symbol %in% lake_dfc_subset_merge_df$gene_symbol)

    donor_ctx_idxs = intersect(which(donorDat$samples$brain == donor), which(donorDat$samples$top_level == 'CTX'))

    ahba_ctx_dat    = 2^donorDat$micro[keep_idxs,donor_ctx_idxs]
    ahba_ctx_probes = donorDat$probes[keep_idxs,]
    ahba_ctx_samples = donorDat$samples[donor_ctx_idxs,]

    ahba_ctx_out = cbind(ahba_ctx_probes$gene_symbol, ahba_ctx_dat)
    colnames(ahba_ctx_out) = c('!Sample_title', ahba_ctx_samples$well_id)

    ciber_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/cibersortX'
    mixture_file = paste0(ciber_dir, '/', region,'_ahba_',donor,'_mixture.txt')
    write(mixture_file,'')
    write.table(x=ahba_ctx_out, file=mixture_file, quote=FALSE, row.names = FALSE, sep='\t')
}








