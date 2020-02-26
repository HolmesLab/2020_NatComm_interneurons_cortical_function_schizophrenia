
library(tidyverse)
library(biomaRt)


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# load snp/gene mapping, created previously
load(verbose=T, paste0(base_dir,'/gene_lists/rr_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.Rdata'))


# read gene ordering pvalb/sst
pvalb_df      = read_csv(paste0(base_dir,'/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
pvalb_ordered = pvalb_df[rev(order(pvalb_df$pvalb_cor)),]

sst_df      = read_csv(paste0(base_dir,'/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))
sst_ordered = sst_df[rev(order(sst_df$sst_cor)),]



# write lists of SNPS tied to the top 500 genes tied to pvalb/sst (for input to GCTA)
lower=1
upper=500

# pvalb SNP list
pvalb_ids  = pvalb_ordered$entrez_id[lower:upper]
pvalb_rsid = as.character(unlist(bim_maf05_gene_snps[as.character(pvalb_ids)]))
pvalb_out  = paste0(base_dir,'/gene_lists/rsid_lists/rr_ahba_rsid_PVALB_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_', as.character(lower), '_', as.character(upper), '.txt')
write.table(file=pvalb_out, x=unique(pvalb_rsid), col.names=F, row.names=F, quote=F)

# sst SNP list
sst_ids    = sst_ordered$entrez_id[lower:upper]
sst_rsid   = as.character(unlist(bim_maf05_gene_snps[as.character(sst_ids)]))
sst_out    = paste0(base_dir,'/gene_lists/rsid_lists/rr_ahba_rsid_SST_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_', as.character(lower), '_', as.character(upper), '.txt')
write.table(file=sst_out, x=unique(sst_rsid), col.names=F, row.names=F, quote=F)



# get information from biomaRt for AHBA genes
grch37        = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="http://grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes    = c("chromosome_name","start_position", "end_position","strand", "ensembl_gene_id", "hgnc_symbol",'entrezgene_id')
gene_info     = getBM(attributes=attributes, filters='entrezgene_id', values=pvalb_ordered$entrez_id, mart=grch37)
nonmatch      = pvalb_ordered[!pvalb_ordered$entrez_id %in% gene_info$entrezgene,]
gene_info_use = gene_info[which(gene_info$chromosome_name %in% as.character(1:22)),]


# LDSC expects ensembl gene ids, which requires some mapping of AHBA gene names using biomaRt

# genes in the GCTA REML gene bin
id_list = list()
id_list[['PVALB']] = pvalb_ids
id_list[['SST']]   = sst_ids

for (gene in c('PVALB','SST')){

  # entrez IDs for top 500 gene bin
  entrez_ids    = id_list[[gene]]

  # pull the corresponding biomaRt gene dictionary (names, positions, ids, etc)
  bin_gene_dict = gene_info_use[which(gene_info_use$entrezgene_id %in% entrez_ids),]

  # check autosomal genes only
  bin_gene_dict = bin_gene_dict[bin_gene_dict$chromosome_name %in% as.character(1:22),]

  # some entrez ids are duplicated. Deduplicate these genes by matching them to AHBA gene/hgnc names
  bin_dups = bin_gene_dict$entrezgene_id[duplicated(bin_gene_dict$entrezgene_id)]
  rm_df    = NULL
  for (dup in bin_dups){
    dup_rows = bin_gene_dict %>% filter(entrezgene_id == dup)
    if (gene == 'PVALB'){
      gene_match = which(!dup_rows$hgnc_symbol %in% pvalb_ordered$gene_symbol[1:500])
    } else if (gene == 'SST'){
      gene_match = which(!dup_rows$hgnc_symbol %in% sst_ordered$gene_symbol[1:500])
    }
    keep_dup   = dup_rows[gene_match,]
    rm_df      = rbind(rm_df, keep_dup)
    print(dup_rows)
  }
  # remove the genes with hgnc names that didnt map to AHBA
  bin_gene_dict = bin_gene_dict %>% filter(!ensembl_gene_id %in% rm_df$ensembl_gene_id)

  # check duplicates
  length(unique(bin_gene_dict$entrezgene_id))
  length(unique(bin_gene_dict$ensembl_gene_id))

  # write ensembl gene lists to be input into LDSC
  bin_ldsc = bin_gene_dict[c('ensembl_gene_id','chromosome_name','start_position','end_position')]
  colnames(bin_ldsc) = c('GENE','CHR','START','END')

  out_file = paste0(base_dir, '/data/ukb/imputed/ldsc_part/GCTA_',as.character(gene),'_ahba_top500_ensembl.txt')
  write.table(bin_ldsc$GENE, out_file, quote=F, row.names=F, col.names=F)
}


# LDSC requires a gene coordinate file to map SNPs to genes in their annotation file format

# format biomaRt table and file with all AHBA genes and chromosomal positions
ldsc_out = gene_info_use[c('ensembl_gene_id', 'chromosome_name', 'start_position', 'end_position')]
colnames(ldsc_out) = c('GENE', 'CHR', 'START', 'END')
ldsc_out = ldsc_out[ldsc_out$CHR %in% as.character(1:22),]
ldsc_out = ldsc_out[order(ldsc_out$CHR),]
ldsc_out = ldsc_out[order(as.numeric(ldsc_out$CHR)),]
head(ldsc_out)

# write gene coordinate file
out_path = paste0(base_dir, '/data/ukb/imputed/ldsc_part/ahba_hg37_gene_coords.txt')
write_delim(ldsc_out, out_path, delim='\t')





# end



