library(tidyverse)
library(biomaRt)
library(EWCE)


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# read bim file containing SNPs analyzed
bim_file      = paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.bim')
bim_maf05_dat = read_delim(bim_file, col_names = c('chr','rsid','V2','bp','a1','a2'), delim='\t')



# Create mart object to get grch37 build snp information
grch37        = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="http://grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes    = c("chromosome_name","start_position", "end_position","strand", "ensembl_gene_id", "hgnc_symbol",'entrezgene_id')
human_ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_dict     = getBM(attributes=attributes, values='*', mart=human_ensembl)



# GTEx lookup table - needed for cross-referencing UKB SNPs to GTEx eQTL snps
gtex_lookup_path = paste0(base_dir, '/data/GTEx/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt')
gtex_lookup      = read_delim(gtex_lookup_path, delim='\t')



# NIH GTEx - read data
gtex_df   = NULL
gtex_regs = c('Anterior_cingulate_cortex_BA24', 'Cortex', 'Frontal_Cortex_BA9')
for (reg in gtex_regs){
  print(reg)
  # read significant pairs
  gtex_path   = paste0(base_dir, '/data/GTEx/GTEx_Analysis_v7_eQTL/Brain_', reg,'.v7.signif_variant_gene_pairs.txt')
  df          = read_delim(gtex_path, delim='\t')
  ens_list    = unlist(lapply(df$gene_id, function(x){strsplit(x,'[.]')[[1]][1] }))
  df$ensemble = ens_list

  # merge snp id with SNP information dataframe
  gtex_lookup_subset = gtex_lookup[gtex_lookup$variant_id %in% df$variant_id,]
  df = merge(df, gtex_lookup_subset, by.x='variant_id', by.y='variant_id')

  # store in dataframe
  gtex_df[[paste0('GTEx_',reg)]] = df
}


# CommonMind - read data
cm_path = paste0(base_dir, '/data/CommonMind/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA-binned.txt')
cm_df   = read_delim(cm_path, delim=' ')
cm_sig  = cm_df[cm_df$FDR %in% c('<0.01', '<0.05'),] # only keep FDR sig snps
cm_sig  = cm_sig[cm_sig$SNP %in% bim_dat$rsid,] # only keep SNPs that are in UKB genotyped data


# the GTEx snp data is large, so it improves lookup speed to split the data by chromosome
bim_maf05_by_chr = NULL
for (chr in 1:22){
  print(chr)
  print(dim(as.data.frame(bim_dat[bim_dat$chr == chr,])))
  bim_maf05_by_chr[[as.character(chr)]] = as.data.frame(bim_dat[bim_dat$chr == chr,])
}

# read gene ordering pvalb/sst
ahba_pvalb = read_csv(paste0(base_dir, '/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
ahba_sst   = read_csv(paste0(base_dir, '/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))

# get information from biomaRt for AHBA genes
gene_info     = getBM(attributes=attributes, filters='entrezgene_id', values=ahba_pvalb$entrez_id, mart=grch37)
nonmatch      = ahba_pvalb[!ahba_pvalb$entrez_id %in% gene_info$entrezgene,]
gene_info_use = gene_info[which(gene_info$chromosome_name %in% as.character(1:22)),]


get_snp_list = function(gene_info, bim_by_chr, bim_dat, gene_snps, ahba_sst, gtex_df, cm_sig){

  # pull intragenic/eQTL snps for each gene
  gene_snps   = list()
  entrez_list = unique(gene_info$entrezgene)

  # retrieve SNP information for each entrez
  ct = 1
  for (entrez in entrez_list ){
    row_idxs = which(gene_info$entrezgene %in% entrez)
    cur_gene_info = gene_info[row_idxs,]

    # if more than one entry is found
    if (nrow(cur_gene_info) > 1){
      print('a')
      use_this_entry = which(cur_gene_info$chromosome_name %in% as.character(1:22))

      # if there are multiple entries still for the entrezid, match based on gene name
      if (length(use_this_entry) > 1){
        entrez_info   = cur_gene_info[use_this_entry,]
        cur_ahba_info = pvalb_df[pvalb_df$entrez_id == entrez,]
        matched_entry = entrez_info[entrez_info$hgnc_symbol %in% cur_ahba_info$gene_symbol,]

        # skip if there is inconsistency between gene symbol ids (hard to be confident of the mapping between AHBA probe and string positions)
        if (nrow(matched_entry) != 1){
          print('skip')
          next
        }
      }
    } else {
      matched_entry = cur_gene_info
    }

    chr       = matched_entry$chromosome_name #gene_info[row,]$chromosome_name
    start_pos = matched_entry$start_position # gene_info[row,]$start_position
    end_pos   = matched_entry$end_position #gene_info[row,]$end_position

    # identify SNPs that fall within start/stop, give or take
    window = 5000
    # if chromosome is NA
    if ( length(chr) == 0 ){
      chr = 'NA'
    }
    # skip if not an autosomal gene
    if ( chr %in% as.character(1:22) == FALSE ){
      next
    }
    cur_df = bim_by_chr[[chr]]

    # snps falling within the intragenic window
    snp_dat = filter(cur_df, bp >= start_pos-window & bp <= end_pos+window)
    
    # skip if there are no intragenic SNPs
    if (nrow(snp_dat) == 0){
      next
    }

    # Get cortical GTEx SNPs
    gtex_rsids = NULL
    for (reg in gtex_regs){
      # pull significant egenes from this regions
      df     = gtex_df[[paste0('GTEx_',reg)]]
      cur_df = df %>% filter(ensemble == matched_entry$ensembl_gene_id)

      # GTEx snps
      cur_rsid   = cur_df$rs_id_dbSNP147_GRCh37p13
      gtex_rsids = c(gtex_rsids, cur_rsid)
    }

    # CommonMind
    cm_top     = cm_sig[which(cm_sig$Gene %in% matched_entry$ensembl_gene_id),]
    gtex_rsids = c(gtex_rsids, cm_top$SNP)

    entrez_name = as.character(entrez)
    entrez_snps = c(snp_dat$rsid, gtex_rsids)
    entrez_snps = entrez_snps[entrez_snps %in% bim_dat$rsid]
    gene_snps[[entrez_name]] = unique(entrez_snps)

    write(paste0(as.character(ct),'/', length(entrez_list)),'')
    ct = ct + 1
  }
  return(gene_snps)
}


# read gene ordering pvalb/sst
pvalb_df = read_csv(paste0(base_dir,'/gene_lists/pvalb_ztransform_ahba_ctx_correlations.csv'))
sst_df   = read_csv(paste0(base_dir,'/gene_lists/sst_ztransform_ahba_ctx_correlations.csv'))



bim_maf05_gene_snps = get_snp_list(gene_info=gene_info,
                                    bim_by_chr=bim_maf05_by_chr,
                                    bim_dat=bim_maf05_dat,
                                    ahba_sst=ahba_sst,
                                    gtex_df=gtex_df,
                                    cm_sig=cm_sig)

# bim_maf05_gene_snps = get_snp_list(bim_by_chr=bim_maf05_by_chr, ahba_sst=ahba_sst, hg37_reannotate=hg37_reannotate, hg37_refgene=hg37_refgene, gtex_df=gtex_df, cm_sig=cm_sig)
save(x=bim_maf05_gene_snps, file=paste0(base_dir,'/gene_lists/rr_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.Rdata'))
load(paste0(base_dir,'/gene_lists/rr_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.05_n10313_rm025.Rdata'), verbose=T)





# read bim file containing SNPs analyzed
bim_maf01_dat = read_delim(paste0(base_dir, '/data/ukb/genotyped/ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.01_n10313.bim'),
                      col_names = c('chr','rsid','V2','bp','a1','a2'), delim='\t')

# the GTEx snp data is large, so it improves lookup speed to split the data by chromosome
bim_maf01_by_chr = NULL
for (chr in 1:22){
  print(chr)
  dim(as.data.frame(bim_maf01_dat[bim_maf01_dat$chr == chr,]))
  bim_maf01_by_chr[[as.character(chr)]] = as.data.frame(bim_maf01_dat[bim_maf01_dat$chr == chr,])
}

bim_maf01_gene_snps = get_snp_list(gene_info=gene_info, bim_by_chr=bim_maf01_by_chr, bim_dat=bim_maf01_dat, ahba_sst=ahba_sst, gtex_df=gtex_df, cm_sig=cm_sig)
save(x=bim_maf01_gene_snps, file=paste0(base_dir,'/gene_lists/rr_SNP_LIST_ukb_geno_v3_caucasian_geno0.02_mind0.1_hwe1e-6_maf0.01_n10313.Rdata'))






