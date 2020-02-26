library(tidyverse)


# Set up directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# read gene ordering pvalb/sst
pvalb_df = read_csv(paste0(base_dir,'/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
sst_df   = read_csv(paste0(base_dir,'/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))


# magma gene-based enrichment
magma_genes = read.table(paste0(base_dir, '/data/scz_gwas/snp_loc_ckqny.scz2snpres.genes.out'), header = T)


# Cross-reference MAGMA/AHBA genes (not all AHBA genes have MAGMA effect estimates)
pvalb_df_match   = pvalb_df[pvalb_df$entrez_id %in% magma_genes$GENE,]
pvalb_magma_sort = pvalb_df_match[rev(order(pvalb_df_match$pvalb_cor)),]

sst_df_match     = sst_df[sst_df$entrez_id %in% magma_genes$GENE,]
sst_magma_sort   = sst_df_match[rev(order(sst_df_match$sst_cor)),]

pvalb_magma_sort$entrezgene = pvalb_magma_sort$entrez_id
sst_magma_sort$entrezgene   = sst_magma_sort$entrez_id


# make sure there is zero overlap between gene lists
intersect(pvalb_magma_sort$entrez_id[1:500], sst_magma_sort$entrez_id[1:500])


# write a gene-list to be used in the partitioned MAGMA analyses (bin size=500)
bottom=1
top=500
out_file = paste0(base_dir, '/data/magma/gene_partitions/ahba_sst_pvalb.txt')
sink(out_file)
while ( top <= 10000 ){
  sst_out   = c(paste0('sst_', as.character(bottom), '_', as.character(top), ' ', paste0(as.character(sst_magma_sort$entrezgene[bottom:top]), collapse=' ')))
  pvalb_out = c(paste0('pvalb_', as.character(bottom), '_', as.character(top), ' ', paste0(as.character(pvalb_magma_sort$entrezgene[bottom:top]), collapse=' ')))
  cat(paste0(sst_out, '\n'))
  cat(paste0(pvalb_out, '\n'))
  bottom=bottom+500
  top=top+500
}
sink()


# Run competitive gene set analysis
cmd = paste( paste0(base_dir,'/external/magma'),
              '--gene-results', paste0(base_dir, '/data/scz_gwas/snp_loc_ckqny.scz2snpres.genes.raw'),
              '--set-annot', paste0(base_dir, '/data/magma/gene_partitions/ahba_sst_pvalb.txt'),
              '--out', paste0(base_dir, '/data/magma/gene_partitions/ahba_sst_pvalb'))
system(cmd)


# read MAGMA results
sst_pvalb_sets = read.table(header=T, paste0(base_dir, '/data/magma/gene_partitions/ahba_sst_pvalb.gsa.out'))


sst_pvalb_sets$SET  = as.character(sst_pvalb_sets$VARIABLE)
sst_pvalb_sets$gene = unlist(lapply(sst_pvalb_sets$SET, function(x) strsplit(x,'_' )[[1]][1]))
sst_pvalb_sets$partition = unlist(lapply(sst_pvalb_sets$SET, function(x) paste0(strsplit(x,'_' )[[1]][2:3], collapse='_')))
sst_pvalb_sets$partition = factor(sst_pvalb_sets$partition, levels=unique(sst_pvalb_sets$partition))


# split into pvalb/sst sets
pvalb_subset = sst_pvalb_sets %>% filter(gene=='pvalb')
cor.test(pvalb_subset$BETA, 1:length(pvalb_subset$BETA))
cor.test(pvalb_subset$BETA, 1:length(pvalb_subset$BETA), method='spearman')
sst_subset = sst_pvalb_sets %>% filter(gene=='sst')
cor.test(sst_subset$BETA, 1:length(sst_subset$BETA), method='spearman')


# plot PVALB rolling window
mypal = rev(colorRampPalette( c( "white", '#EBF0B5', '#46B6C2', '#283B8E', '#283B8E') )( 20 ))
pvalb_subset$idx  = 1:length(pvalb_subset$BETA)
pvalb_subset$line = predict(lm('BETA ~ idx', data=pvalb_subset))

pdf(paste0(base_dir, '/figures/PaperFig_pvalb_magma_enrichment.pdf'), height=2, width=3.5)
ggplot(pvalb_subset, aes(idx, BETA)) +
  geom_bar(aes(fill = partition), stat = "identity", position = "dodge") +
  geom_line(aes(idx, line), linetype = 2) +
  scale_y_continuous(breaks = seq(-15,.15,.05), expand = c(0, 0), limits=c(-.15,.15)) +
  scale_fill_manual(values = mypal ) +
  theme_classic()
dev.off()


# plot SST rolling window
mypal = rev(colorRampPalette( c( "white", '#FCEDA9', '#EF924F', '#AE232D', '#AE232D') )( 20 ))
sst_subset$idx  = 1:length(sst_subset$BETA)
sst_subset$line = predict(lm('BETA ~ idx', data=sst_subset))

pdf(paste0(base_dir, '/figures/PaperFig_sst_magma_enrichment.pdf'), height=2, width=3.5)
ggplot(sst_subset, aes(idx, BETA)) +
  geom_bar(aes(fill = partition), stat = "identity", position = "dodge") +
  geom_line(aes(idx, line), linetype = 2) +
  scale_y_continuous(breaks = seq(-15,.15,.05), expand = c(0, 0), limits=c(-.15,.15)) +
  scale_fill_manual(values = mypal ) +
  theme_classic()
dev.off()




# plot spectrum for making color bars
for (gene in c('sst','pvalb')){
  if (gene == 'sst'){
    mypal <- colorRampPalette( c( "white", "white", '#FCEDA9', '#EF924F', '#AE232D', '#AE232D') )( 100 )
  } else if (gene == 'pvalb'){
    mypal <- colorRampPalette( c( "white", "white", '#EBF0B5', '#46B6C2', '#283B8E', '#283B8E') )( 100 )
  }
  # color scale for figure 1
  z=matrix(1:100,nrow=1)
  x=1
  y=seq(-3,3,len=100)
  image(x,y,z,col=mypal,axes=FALSE,xlab="",ylab="")
}










