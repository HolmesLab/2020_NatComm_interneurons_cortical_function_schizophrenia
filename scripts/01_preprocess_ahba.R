# Load packages
library(tidyverse)
library(WGCNA)


# this function will map AHBA samples to their overarching regional category (e.g. mPFC >> CORTEX; mediodorsal_thal >> THAL)
find_top_level = function(ref_df, struct_id, ontology, top_level) {

    ontology_row  = which(ontology$id == struct_id )

    # print feedback if more than one match, this should never be the case
    if (length(ontology_row) > 1){
        write("ERROR")
    }

    ontology_info = ontology[ontology_row,]
    splits    = strsplit(as.character(ontology_info$structure_id_path), '/')[[1]]
    reg.match = ref_df$id[which(ref_df$id %in% splits)]
    if (length(reg.match) == 1){
        region_cat   = as.character(ref_df$top_level[which(ref_df$id %in% splits)])
        region_name  = ref_df[ref_df$id == reg.match,]$name
        return(c(reg.match, region_cat, region_name))
    } else {
        return(c(NA,NA,NA))
    }
}


# Modify these filepaths for your local directory structure
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'


# AHBA subject list/path to AHBA microarray data
data_path  = paste(base_dir, '/data/ahba', sep = '') 
donor_nums = c('9861', '10021', '12876', '14380', '15496', '15697')


# Map ontology structure names to more general sample categories
# e.g. left motor gyrus --> CTX
top_level = read_csv(paste0(base_dir, '/ref_files/top_level_categories.csv'))


# Read data from each subject
ahba_data  = NULL
sampleInfo = NULL
microExprA = matrix()
for ( donor in donor_nums ) {
    file = paste('donor', donor, sep='')
    write(paste('Reading and collapsing data for donor: ', donor, sep = ''),'')

    # create donor field
    ahba_data[[donor]] = NULL


    # read sample Information
    saname    = paste(data_path, file, 'SampleAnnot.csv', sep='/')
    samp_info = read_csv(saname)
    ahba_data[[donor]]$raw_samp = samp_info
    samp_info$brain = donor


    # concatenate sample info into one overall matrix
    sampleInfo = rbind(sampleInfo, samp_info)

    # Gene Ontology Info
    oname    = paste(data_path, file, 'Ontology.csv', sep='/')
    ont_data = read_csv(oname)
    ahba_data[[donor]]$raw_ont = ont_data
    ontology = ont_data

    # identify samples in the regions we want to analyze
    out         = do.call(rbind, lapply(samp_info$structure_id, find_top_level, ref_df=top_level, ontology=ontology))
    reg_info    = as.data.frame(out)
    reg_info$V1 = as.numeric(as.character(reg_info$V1))
    colnames(reg_info) = c('reg_num', 'top_level', 'region_clean')

    # Read Microarray Data
    fname      = paste(data_path, file, 'MicroarrayExpression.csv', sep='/')
    microdata  = read_csv(fname, col_names = F)
    micro_arr  = microdata

    # remove first column contains probe IDs
    micro_temp           = as.data.frame(microdata[,-which(colnames(microdata) == 'X1')]) # First column contains probe IDs
    rownames(micro_temp) = micro_arr$X1 # collapseRows requires rownames
    micro_df             = as.data.frame(micro_temp)
    ahba_data[[donor]]$raw_micros = micro_df

    # create an overall dataframe with all expression data
    microExprA = cbind(microExprA, micro_df)

    # Information about each Gene Probe (~50,000 probes for ~20,000 genes)
    pname     = paste(data_path, file, 'Probes.csv', sep='/')
    probes    = read_csv(pname)
    ahba_data[[donor]]$raw_probes = probes
    probeInfo = probes

    # discard probes without an entrez-id
    num_samples   = dim(ahba_data[[donor]]$raw_micros)[2]
    trash_me      = is.na(ahba_data[[donor]]$raw_probes$entrez_id) 
    good_probes   = ahba_data[[donor]]$raw_probes[trash_me == FALSE,]

    # PA
    write('PA','')
    pname    = paste(data_path, file, 'PACall.csv', sep='/')
    pa_data  = read.csv(pname, header=F)
    rownames(pa_data) = pa_data$V1
    pa_data$V1 = NULL
    ahba_data[[donor]]$raw_pas = pa_data


    # Select the probes with valid Entrez IDs
    ahba_data[[donor]]$probes_filter = good_probes
    ahba_data[[donor]]$micro_filter  = ahba_data[[donor]]$raw_micros[rownames(ahba_data[[donor]]$raw_micros) %in% good_probes$probe_id,]
    ahba_data[[donor]]$pas_filter    = ahba_data[[donor]]$raw_pas[rownames(ahba_data[[donor]]$raw_pas) %in% good_probes$probe_id,]


    ctx_samples = which(reg_info$top_level == 'CTX')
    n_samples   = length(ctx_samples)

    valid_entrez_idxs = which(!is.na(ahba_data[[donor]]$raw_probes$entrez_id))
    ahba_data[[donor]]$raw_ctx_pas = ahba_data[[donor]]$pas_filter[,ctx_samples]
    ahba_data[[donor]]$probe_noise = rowSums(ahba_data[[donor]]$raw_ctx_pas)/n_samples
}
microExprA$microExprA = NULL

# save various forms of the data
# ----

# microarray expression; 58692 probes x 3702 samples
save(x=microExprA, file=paste0(base_dir, '/data/ahba/ahba_microExprA.Rdata'))

# probe metadata; 58692 probes x 7 info columns
save(x=probeInfo, file=paste0(base_dir, '/data/ahba/ahba_probeInfo.Rdata')) 

# sample information; 3702 samples x 14 info columns
save(x=sampleInfo, file=paste0(base_dir, '/data/ahba/ahba_sampleInfo.Rdata'))

# larger data object, with data split by subject
save(x=ahba_data, file=paste0(base_dir, '/data/ahba/ahba_data_object.Rdata'))
load(paste0(base_dir, '/data/ahba/ahba_data_object.Rdata'), verbose=T)



# remove probes without an entrez id

# 48171 x 3702
useMicro  = microExprA[!is.na(probeInfo$entrez_id),]
# 48171 x 7
useProbes = probeInfo[!is.na(probeInfo$entrez_id),]


# average probe noise across subjects. 
# probe noise is number of samples with significant non-zero expression, divided by total samples 
avg_probe_noise = rowMeans(cbind(ahba_data[[1]]$probe_noise,
                            ahba_data[[2]]$probe_noise,
                            ahba_data[[3]]$probe_noise,
                            ahba_data[[4]]$probe_noise,
                            ahba_data[[5]]$probe_noise,
                            ahba_data[[6]]$probe_noise))


# keep probes that are, on average across subs, expressed in at least 20% of cortical samples
useMicro  = useMicro[which(avg_probe_noise > .2),] # 36762 x 3702
useProbes = useProbes[which(avg_probe_noise > .2),] # 36762 x 7


# map each ahba subregion to its overall category
out         = do.call(rbind, lapply(sampleInfo$structure_id, find_top_level, ref_df=top_level, ontology=ontology))
reg_info    = as.data.frame(out)
reg_info$V1 = as.numeric(as.character(reg_info$V1))
colnames(reg_info) = c('reg_num', 'top_level', 'region_clean')


# add new sample super-category info to original sample information dataframe
sampleInfo    = cbind(sampleInfo, reg_info)
micro_subset  = useMicro[,which(!is.na(sampleInfo$top_level))]
sample_subset = sampleInfo[which(!is.na(sampleInfo$top_level)),]

write_csv(sampleInfo, path=paste0(base_dir, '/data/ahba/ahba_sampleInfo.csv'))


# Collapse probes based on highest mean expression
rownames(micro_subset) = as.character(useProbes$probe_name)
datExpr           = collapseRows(datET=micro_subset, rowGroup=as.character(useProbes$gene_symbol), rowID=as.character(useProbes$probe_name), method="MaxMean", connectivityBasedCollapsing=FALSE)
donorDat = NULL
donorDat$micro    = datExpr$datETcollapsed
donorDat$samples  = sample_subset
donorDat$ontology = ontology
probes2keep       = as.character(datExpr$group2row[,2])
probesKept        = probeInfo[probeInfo$probe_name %in% probes2keep,]
donorDat$probes   = probesKept[order(probesKept$gene_symbol),]
length(which(donorDat$probes$probe_name == probes2keep))


# For quicker loading after the above preprocessing has been run, save Rdata structure
save(x=donorDat, file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'))


