library(tidyverse)
library(gplots)
library(RColorBrewer)
library(Cairo)
library(cifti)
library(cowplot)


# function to plot parcel-wise cortical correlations for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){
    base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
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
base_dir  = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
ciber_dir = paste0(base_dir, '/data/cibersortX')


# read cell type expression of "signature genes"
vis_gep_file = paste0(ciber_dir, '/CIBERSORTx_Job1_Lake_VisualCortex_ahba_matched_sc_sample_inferred_phenoclasses.CIBERSORTx_Job1_Lake_VisualCortex_ahba_matched_sc_sample_inferred_refsample.bm.K999.txt')
vis_gep = read.table(vis_gep_file, header=T)


# visual cortex - create heatmap of gene signatures
file_out = paste0(paste0(base_dir, '/figures/vis_gep_signature_matrix.pdf'))
CairoPDF(file_out, height=5, width=4)
heatmap.2(as.matrix(vis_gep[,2:ncol(vis_gep)]), col=brewer.pal(11,"RdBu"), key=FALSE, dendrogram='col', scale="row", trace="none", labRow=FALSE)
dev.off()
file_out




# average expression of signature genes in each cell class
vis_sig_cormat = cor(vis_gep[,2:ncol(vis_gep)], use = 'pairwise.complete')
cell_text = t(as.matrix(vis_sig_cormat))
cell_text = round(cell_text,2)
# color pallette for the correlation plot
my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
my_palette = rev(my_palette)

# plot
pdf(file = paste0(base_dir, '/figures/lake_vis_corr_of_signature_mat.pdf'), width=10, height=10)
heatmap.2(as.matrix(vis_sig_cormat), cellnote=cell_text, trace="none", notecol='black',col=my_palette,breaks=seq(-1, 1, length.out=300),
            distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
            hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm=F,symkey=F,symbreaks=T)
dev.off()


# signature gene expression data
dfc_gep_file = paste0(ciber_dir, '/CIBERSORTx_Job1_Lake_FrontalCortex_ahba_matched_sc_sample_inferred_phenoclasses.CIBERSORTx_Job1_Lake_FrontalCortex_ahba_matched_sc_sample_inferred_refsample.bm.K999.txt')
dfc_gep      = read.table(dfc_gep_file, header=T)


# frontal cortex - create heatmap of gene signatures
file_out = paste0(paste0(base_dir, '/figures/PaperFig_dfc_gep_signature_matrix.pdf'))
CairoPDF(file_out, height=5, width=4)
heatmap.2(as.matrix(dfc_gep[,2:ncol(dfc_gep)]), col=brewer.pal(11,"RdBu"), key=FALSE, dendrogram='col', scale="row", trace="none", labRow=FALSE)
dev.off()
file_out


# average expression of signature genes in each cell class
dfc_sig_cormat = cor(dfc_gep[,2:ncol(dfc_gep)], use = 'pairwise.complete')
cell_text = t(as.matrix(dfc_sig_cormat))
cell_text = round(cell_text,2)
pdf(file = paste0(base_dir, '/figures/lake_dfc_corr_of_signature_mat.pdf'), width=10, height=10)
heatmap.2(as.matrix(dfc_sig_cormat), cellnote=cell_text, trace="none", notecol='black',col=my_palette,breaks=seq(-1, 1, length.out=300),
            distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
            hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm=F,symkey=F,symbreaks=T)
dev.off()



# read AHBA information
load(file=paste0(base_dir, '/data/ahba/ahba_data_object.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)
donor_arr = names(ahba_data)

# read sample-to-vertex projection info
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))
sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]

donor_arr = c("9861","10021","12876","14380","15496","15697")

# read cibersort data from each donor and combine
for ( region in c('FrontalCortex','VisualCortex') ){

    write(region,'')
    all_ciber = NULL
    for (donor in donor_arr){
        write(donor, '')

        path_to_ciber = paste0(ciber_dir, '/',region,'_ahba_',donor,'_Fractions.txt')
        ciber_dat = read_delim(path_to_ciber, delim='\t') #read.table(path_to_ciber, header=T)
        all_ciber = rbind(all_ciber, ciber_dat)

        ciber_field = paste0(region, '_CiberFrac')
        ahba_data[[donor]][[ciber_field]] = merge(x=ciber_dat, y=sample_info, by.y='well_id', by.x='Mixture')
    }
    write_csv(x=all_ciber, path=paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_',region,'.csv'))
}


# Read cell fractions derived from Visual Cortex cell signatures
vis_ciber = read_csv(paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_VisualCortex.csv'))
colnames(vis_ciber) = paste0('VIS_',colnames(vis_ciber))


# Read cell fractions derived from Frontal Cortex cell signatures
dfc_ciber = read_csv(paste0(ciber_dir, '/Cibersortx_AHBAsamples_Lake_FrontalCortex.csv'))
colnames(dfc_ciber) = paste0('DFC_', colnames(dfc_ciber))


# frontal cortex spatial correlations
dfc_cormat   = as.data.frame(cor(dfc_ciber[colnames(dfc_ciber)[2:19]], use='pairwise.complete'))
dfc_sst_df   = data.frame(sst_cors = dfc_cormat$DFC_SST, cell = rownames(dfc_cormat), region='DFC')
dfc_pvalb_df = data.frame(pvalb_cors = dfc_cormat$DFC_PVALB, cell = rownames(dfc_cormat), region='DFC')


# visual cortex spatial correlations
vis_cormat = as.data.frame(cor(vis_ciber[colnames(vis_ciber)[2:19]], use='pairwise.complete'))
vis_sst_df = data.frame(sst_cors = vis_cormat$VIS_SST, cell = rownames(vis_cormat), region='VIS')
vis_pvalb_df = data.frame(pvalb_cors = vis_cormat$VIS_PVALB, cell = rownames(vis_cormat), region='VIS')


# Plot Frontal Cortex spatial correlations to SST
celltypecors = dfc_sst_df[order(dfc_sst_df$sst_cors),]
celltypecors = celltypecors[!is.na(celltypecors$sst_cors),]
celltypecors = celltypecors[celltypecors$sst_cors != 1,]
celltypecors$cell = gsub('DFC_', '', celltypecors$cell)
celltypecors$cell = factor(as.character(celltypecors$cell), levels=as.character(celltypecors$cell))

# Frontal Cortex spatial correlations to SST cell fraction
gg_dfc_sst = ggplot(data=celltypecors, aes(x=cell, y=sst_cors)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        #scale_fill_manual(values=color_arr[color_classes %in% celltypecors$cell_color]) +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-.6, .6, by=.2),2), limits=c(-.6,.6))


# Frontal Cortex spatial correlations to PVALB cell fraction
celltypecors = dfc_pvalb_df[order(dfc_pvalb_df$pvalb_cors),]
celltypecors = celltypecors[!is.na(celltypecors$pvalb_cors),]
celltypecors = celltypecors[celltypecors$pvalb_cors != 1,]
celltypecors$cell = gsub('DFC_', '', celltypecors$cell)
celltypecors$cell = factor(as.character(celltypecors$cell), levels=as.character(celltypecors$cell))

gg_dfc_pvalb = ggplot(data=celltypecors, aes(x=cell, y=pvalb_cors)) +
    geom_bar(stat="identity") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-.6, .6, by=.2),2), limits=c(-.6,.6))



# Plot Visual Cortex spatial correlations to PVALB
celltypecors = vis_pvalb_df[order(vis_pvalb_df$pvalb_cors),]
celltypecors = celltypecors[!is.na(celltypecors$pvalb_cors),]
celltypecors = celltypecors[celltypecors$pvalb_cors != 1,]
celltypecors$cell = gsub('VIS_', '', celltypecors$cell)
celltypecors$cell = factor(as.character(celltypecors$cell), levels=as.character(celltypecors$cell))

gg_vis_pvalb = ggplot(data=celltypecors, aes(x=cell, y=pvalb_cors)) +
    geom_bar(stat="identity") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-.6, .6, by=.2),2), limits=c(-.6,.6))



# Plot Visual Cortex spatial correlations to SST
celltypecors = vis_sst_df[order(vis_sst_df$sst_cors),]
celltypecors = celltypecors[!is.na(celltypecors$sst_cors),]
celltypecors = celltypecors[celltypecors$sst_cors != 1,]
celltypecors$cell = gsub('VIS_', '', celltypecors$cell)
celltypecors$cell = factor(as.character(celltypecors$cell), levels=as.character(celltypecors$cell))

gg_vis_sst = ggplot(data=celltypecors, aes(x=cell, y=sst_cors)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-.6, .6, by=.2),2), limits=c(-.6,.6))


figure_out = paste0(base_dir,'/figures/PaperFig_celltype_VisDfc_cor_to_sst.pdf')
CairoPDF(figure_out, width=3, height=3)
plot_grid(gg_dfc_sst, gg_vis_sst)
dev.off()

figure_out = paste0(base_dir,'/figures/PaperFig_celltype_VisDfc_cor_to_pvalb.pdf')
CairoPDF(figure_out, width=3, height=3)
plot_grid(gg_dfc_pvalb, gg_vis_pvalb)
dev.off()


# merge cibersort data with AHBA sample dat
cell_types = colnames(dfc_ciber)[2:19]
cell_types = cell_types[!cell_types %in% c('DFC_Ex6')]
ciber_comb = merge(x=dfc_ciber, y=vis_ciber, by.x='DFC_Mixture', by.y='VIS_Mixture')
sample_dat_ciber = merge(x=sample_dat, y=ciber_comb, by.x='well_id', by.y='DFC_Mixture', all.x=T)


dim(sample_dat_ciber)
head(sample_dat_ciber)


# Collapse individual cell proportions into schaeffer parcels

# stitch the left/right verticies together to match Schaeffer parcel cifti format
sample_dat_ciber$bihemi_vertex = sample_dat_ciber$vertex + 1 # cifti indices index at 0, R indexes at 1
right_ctx_idxs = intersect(grep('right', sample_dat_ciber$structure_name), which(sample_dat_ciber$top_level == 'CTX'))
sample_dat_ciber$bihemi_vertex[right_ctx_idxs] = sample_dat_ciber$bihemi_vertex[right_ctx_idxs] + 32492


parcel_average = function(ciber_dat, ref_region, parcels, net, cell_types, name_string=''){
    donor_specific_expression = NULL

    # schaeffer parcellation by vertex
    schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
    schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)

    # corresponding parcel labels
    schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
    schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]

    # calculate the average gene expression within each parcel
    sst_pvalb_arr = NULL
    schaeffer_mat = matrix(NA, ncol=as.numeric(parcels), nrow=length(cell_types))
    for (donor in donor_arr){donor_specific_expression[[donor]] = schaeffer_mat}

    reg_cell_types = paste0(ref_region, '_', cell_types)
    #reg_ciber = ciber_dat[reg_cell_types]

    for (idx in 1:length(schaef_labels)){
        write(idx,'')
        parcel_idxs   = which(schaef_cii$data == idx) # schaeffer indices
        match_idxs    = which(ciber_dat$bihemi_vertex %in% parcel_idxs) # foci within this parcel
        match_samples = ciber_dat[match_idxs,] # data for this parcel

        schaeffer_expr = colMeans(ciber_dat[match_idxs, reg_cell_types], na.rm=T) # average expressino of every gene, across samples in this parcel
        schaeffer_mat[,idx] = schaeffer_expr # plug in values to the pre-allocated matrix

        for (donor in donor_arr){
            donor_idxs     = which(as.character(ciber_dat$brain) == donor)
            donor_reg_idxs = intersect(match_idxs, donor_idxs)
            expr = colMeans(ciber_dat[donor_reg_idxs,reg_cell_types]) # average expressino of every gene, across samples in this parcel
            donor_specific_expression[[donor]][,idx] = expr
        }
    }
    schaef_out = as.data.frame(schaeffer_mat)
    rownames(schaef_out) = reg_cell_types
    colnames(schaef_out) = schaef_labels

    for (donor in donor_arr){
        donor_specific_expression[[donor]] = as.data.frame(donor_specific_expression[[donor]])
        rownames(donor_specific_expression[[donor]]) = reg_cell_types
        colnames(donor_specific_expression[[donor]]) = schaef_labels

        donor_specific_expression[[donor]]$gene = reg_cell_types
        write_csv(donor_specific_expression[[donor]], path=paste0(base_dir, '/data/ahba/schaeffer_LAKE_',ref_region,'_donor',as.character(donor),'_p', parcels,'_',net,'Net_expr_mat',name_string,'.csv'))
    }

    schaef_out$gene = reg_cell_types
    write_csv(schaef_out, path=paste0(base_dir, '/data/ahba/schaeffer_LAKE_',ref_region,'_',parcels,'_',net,'Net_expr_mat',name_string,'.csv'))

    return(list(schaef_out, donor_specific_expression))
}

# average parcel-wise fraction of cell types (DFC)
cell_types = gsub('DFC_', '', cell_types)
cell_types = colnames(dfc_ciber)[2:19]
cell_types = cell_types[!cell_types %in% c('DFC_Ex6')]
cell_types = gsub('DFC_', '', cell_types)

# do it
DFC_400_7_schaef_out = parcel_average(ciber_dat=sample_dat_ciber, ref_region='DFC', parcels='400', net='7', cell_types=cell_types)
DFC_400_7_donor_specific_expression = DFC_400_7_schaef_out[[2]]
DFC_400_7_schaef_out = DFC_400_7_schaef_out[[1]]


# average parcel-wise fraction of cell types (VIS)
cell_types = colnames(vis_ciber)[2:19]
cell_types = cell_types[!cell_types %in% c('VIS_Ex2','VIS_Ex6')]
cell_types = gsub('VIS_', '', cell_types)

# do it
VIS_400_7_schaef_out = parcel_average(ciber_dat=sample_dat_ciber, ref_region='VIS', parcels='400', net='7', cell_types=cell_types)
VIS_400_7_donor_specific_expression = VIS_400_7_schaef_out[[2]]
VIS_400_7_schaef_out = VIS_400_7_schaef_out[[1]]


# correlate PVALB cell fraction distributions derived from DFC/VIS cell signatures
x = as.numeric(VIS_400_7_schaef_out[VIS_400_7_schaef_out$gene == 'VIS_PVALB',1:400])
y = as.numeric(DFC_400_7_schaef_out[DFC_400_7_schaef_out$gene == 'DFC_PVALB',1:400])
cor.test(x,y)

# correlate SST cell fraction distributions derived from DFC/VIS cell signatures
x = as.numeric(VIS_400_7_schaef_out[VIS_400_7_schaef_out$gene == 'VIS_SST',1:400])
y = as.numeric(DFC_400_7_schaef_out[DFC_400_7_schaef_out$gene == 'DFC_SST',1:400])
cor.test(x,y)





# DFC - create averaged plots
nparcel = '400'
net_num = '7'
lake_dfc_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_DFC_',nparcel,'_',net_num,'Net_expr_mat.csv'))
for (cell in lake_dfc_dat$gene){
    write(cell,'')
    cell_proportion = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == cell,1:as.numeric(nparcel)])
    cell_proportion[which(cell_proportion == 0)] = .001
    plot_matlab(values=as.numeric(cell_proportion), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_Lake_DFC_',cell,'_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)
}
sst_proportion   = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == 'DFC_SST',1:as.numeric(nparcel)])
pvalb_proportion = as.numeric(lake_dfc_dat[lake_dfc_dat$gene == 'DFC_PVALB',1:as.numeric(nparcel)])
sst_minus_pvalb  = scale(sst_proportion) - scale(pvalb_proportion)
sst_minus_pvalb[which(sst_minus_pvalb == 0)] = .001
plot_matlab(values=as.numeric(sst_minus_pvalb), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_schaeffer_LAKE_DFC_SST-PVALB_schaef_',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)



# VIS - create averaged plots
nparcel = '400'
net_num = '7'
lake_vis_dat = read_csv(paste0(base_dir, '/data/ahba/schaeffer_LAKE_VIS_',nparcel,'_',net_num,'Net_expr_mat.csv'))
for (cell in lake_vis_dat$gene){
    write(cell,'')
    cell_proportion = as.numeric(lake_vis_dat[lake_vis_dat$gene == cell, 1:as.numeric(nparcel)])
    cell_proportion[which(cell_proportion == 0)] = .001
    plot_matlab(values=as.numeric(cell_proportion), out_path=paste0(base_dir, '/figures/surface_plots/PaperPlot_Lake_',cell,'_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)
}
sst_proportion   = as.numeric(lake_vis_dat[lake_vis_dat$gene == 'VIS_SST',1:as.numeric(nparcel)])
pvalb_proportion = as.numeric(lake_vis_dat[lake_vis_dat$gene == 'VIS_PVALB',1:as.numeric(nparcel)])
sst_minus_pvalb  = scale(sst_proportion) - scale(pvalb_proportion)
sst_minus_pvalb[which(sst_minus_pvalb == 0)] = .001
plot_matlab(values=as.numeric(sst_minus_pvalb), out_path=paste0(base_dir, '/figures/surface_plots/PaperFig_schaeffer_LAKE_VIS_SST-PVALB_schaef_',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num=nparcel, net_num=net_num)


# Spatial correlation of DFC cell types
dfc_sc_spatial_cors = cor(t(lake_dfc_dat[,1:400]), use = 'pairwise.complete')
colnames(dfc_sc_spatial_cors) = lake_dfc_dat$gene
rownames(dfc_sc_spatial_cors) = lake_dfc_dat$gene


# get ready for plotting
my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
my_palette = rev(my_palette)
cell_text  = t(as.matrix(dfc_sc_spatial_cors))
cell_text  = round(cell_text,2)


# plot corrmat!
pdf(file = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperPlot_lake_dfc_spatial_cors.pdf', width=10, height=10)
heatmap.2(as.matrix(dfc_sc_spatial_cors), cellnote=cell_text, trace="none", notecol='black',col=my_palette,
            distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
            hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm = TRUE)
dev.off()



# Spatial correlation of VIS cell types
x = lake_vis_dat[lake_vis_dat$gene != 'VIS_In3',]
vis_sc_spatial_cors = cor(t(x[,1:400]), use = 'pairwise.complete')
colnames(vis_sc_spatial_cors) = x$gene
rownames(vis_sc_spatial_cors) = x$gene

# get ready for plotting
my_palette = colorRampPalette(c('#67001F','#B2182B','#D6604D','#F4A582','#FDDBC7','#F7F7F7','#D1E5F0','#92C5DE','#4393C3','#2166AC','#053061'))(n = 299)
my_palette = rev(my_palette)
cell_text  = t(as.matrix(vis_sc_spatial_cors))
cell_text  = round(cell_text,2)

# plot corrmat!
pdf(file = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperPlot_lake_vis_spatial_cors.pdf', width=10, height=10)
heatmap.2(as.matrix(vis_sc_spatial_cors), cellnote=cell_text, trace="none", notecol='black',col=my_palette,
            distfun=function(x) dist(x, method="euclidean"),dendrogram="col",density.info="none",
            hclustfun=function(x) hclust(x, method="ward.D2"),Rowv=T, symm = TRUE)
dev.off()




# read parcelwise expression of each individual AHBA gene
nparcel = '400'
net_num = '7'
ahba_expr = read_csv(paste0(base_dir, '/data/ahba/schaeffer',nparcel,'_',net_num,'Net_expr_mat.csv'))
sst_expr   = t(ahba_expr[ahba_expr$gene == 'SST',1:400])
pvalb_expr = t(ahba_expr[ahba_expr$gene == 'PVALB',1:400])


# prep VIS data
lake_vis_just_expr = as.data.frame(t(lake_vis_dat[,1:400]))
lake_vis_just_expr$region = rownames(lake_vis_just_expr)
colnames(lake_vis_just_expr)[1:17] = lake_vis_dat$gene
length(which(rownames(pvalb_expr) == rownames(lake_vis_just_expr)))


# VIS to AHBA SST
vis_cor_to_sstGene = cor(lake_vis_just_expr[1:400,1:17], sst_expr[1:400], use='pairwise.complete')
vis_cor_to_sstGene = data.frame(cor2sst=vis_cor_to_sstGene[!is.na(vis_cor_to_sstGene)], cell=rownames(vis_cor_to_sstGene)[!is.na(vis_cor_to_sstGene)])
vis_cor_to_sstGene = vis_cor_to_sstGene[rev(order(vis_cor_to_sstGene$cor2sst)),]
vis_cor_to_sstGene$cell = factor(gsub('VIS_','',vis_cor_to_sstGene$cell), levels= gsub('VIS_','',vis_cor_to_sstGene$cell))
vis_cor_to_sstGene$region = 'VIS'

# VIS to AHBA PVALB
vis_cor_to_pvalbGene = cor(lake_vis_just_expr[1:400,1:17], pvalb_expr[1:400], use='pairwise.complete')
vis_cor_to_pvalbGene = data.frame(cor2pvalb=vis_cor_to_pvalbGene[!is.na(vis_cor_to_pvalbGene)], cell=rownames(vis_cor_to_pvalbGene)[!is.na(vis_cor_to_pvalbGene)])
vis_cor_to_pvalbGene = vis_cor_to_pvalbGene[rev(order(vis_cor_to_pvalbGene$cor2pvalb)),]
vis_cor_to_pvalbGene$cell = factor(gsub('VIS_','',vis_cor_to_pvalbGene$cell), levels= gsub('VIS_','',vis_cor_to_pvalbGene$cell))
vis_cor_to_pvalbGene$region = 'VIS'



# prep DFC data
lake_dfc_just_expr = as.data.frame(t(lake_dfc_dat[,1:400]))
lake_dfc_just_expr$region = rownames(lake_dfc_just_expr)
colnames(lake_dfc_just_expr)[1:17] = lake_dfc_dat$gene

# DFC to AHBA SST
dfc_cor_to_sstGene = cor(lake_dfc_just_expr[1:400,1:17], sst_expr[1:400], use='pairwise.complete')
dfc_cor_to_sstGene = data.frame(cor2sst=dfc_cor_to_sstGene[!is.na(dfc_cor_to_sstGene)], cell=rownames(dfc_cor_to_sstGene)[!is.na(dfc_cor_to_sstGene)])
dfc_cor_to_sstGene = dfc_cor_to_sstGene[rev(order(dfc_cor_to_sstGene$cor2sst)),]
dfc_cor_to_sstGene$cell = factor(gsub('DFC_','',dfc_cor_to_sstGene$cell), levels= gsub('DFC_','',dfc_cor_to_sstGene$cell))
dfc_cor_to_sstGene$region = 'DFC'

# DFC to AHBA PVALB
dfc_cor_to_pvalbGene = cor(lake_dfc_just_expr[1:400,1:17], pvalb_expr[1:400], use='pairwise.complete')
dfc_cor_to_pvalbGene = data.frame(cor2pvalb=dfc_cor_to_pvalbGene[!is.na(dfc_cor_to_pvalbGene)], cell=rownames(dfc_cor_to_pvalbGene)[!is.na(dfc_cor_to_pvalbGene)])
dfc_cor_to_pvalbGene = dfc_cor_to_pvalbGene[rev(order(dfc_cor_to_pvalbGene$cor2pvalb)),]
dfc_cor_to_pvalbGene$cell = factor(gsub('DFC_','',dfc_cor_to_pvalbGene$cell), levels= gsub('DFC_','',dfc_cor_to_pvalbGene$cell))
dfc_cor_to_pvalbGene$region = 'DFC'


# PVALB - combine all of the correlations for plotting
both_cor_to_pvalbGene = rbind(dfc_cor_to_pvalbGene, vis_cor_to_pvalbGene)
both_cor_to_pvalbGene = rbind(both_cor_to_pvalbGene, data.frame(cor2pvalb=NA, cell='Ex2', region='VIS'))
both_cor_to_pvalbGene = rbind(both_cor_to_pvalbGene, data.frame(cor2pvalb=NA, cell='In2', region='DFC'))
both_cor_to_pvalbGene = rbind(both_cor_to_pvalbGene, data.frame(cor2pvalb=NA, cell='In3', region='VIS'))


# order cells by their average correlation to PVALB
mean_cor = both_cor_to_pvalbGene %>% group_by(cell) %>% summarise(x=mean(cor2pvalb, na.rm=T))
cell_ordering = as.character(mean_cor$cell[order(mean_cor$x)])
both_cor_to_pvalbGene$cell = factor(as.character(both_cor_to_pvalbGene$cell), levels=cell_ordering)

cor2pvalb = ggplot(data=both_cor_to_pvalbGene, aes(x=cell, y=cor2pvalb, fill=region)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1, 1, by = .5)) +
    scale_fill_manual(values=c('#6E727F','#BEC1CC')) +
            theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=8, angle = 90, hjust = 1),
          axis.text.y = element_text(colour="black", size=8),
          axis.title.x = element_text(colour="black", size=8),
          axis.title.y = element_blank(),
            legend.position = "none") +
    coord_flip()


# SST - combine all of the correlations for plotting
both_cor_to_sstGene = rbind(dfc_cor_to_sstGene, vis_cor_to_sstGene)
both_cor_to_sstGene = rbind(both_cor_to_sstGene, data.frame(cor2sst=NA, cell='Ex2', region='VIS'))
both_cor_to_sstGene = rbind(both_cor_to_sstGene, data.frame(cor2sst=NA, cell='In2', region='DFC'))
both_cor_to_sstGene = rbind(both_cor_to_sstGene, data.frame(cor2sst=NA, cell='In3', region='VIS'))


# order cells by their average correlation to SST
mean_cor = both_cor_to_sstGene %>% group_by(cell) %>% summarise(x=mean(cor2sst, na.rm=T))
cell_ordering = as.character(mean_cor$cell[rev(order(mean_cor$x))])
both_cor_to_sstGene$cell = factor(as.character(both_cor_to_sstGene$cell), levels=cell_ordering)

cor2sst = ggplot(data=both_cor_to_sstGene, aes(x=cell, y=cor2sst, fill=region)) +
    geom_col(colour="black",width=0.7,
           position=position_dodge(0.7)) +
    #geom_bar(stat="identity", color="black", position=position_dodge(), bar_spacing=1) +
    #geom_errorbar(aes(ymin=norm_expr-expr_se, ymax=norm_expr+expr_se), width=.4, position=position_dodge(.9)) +
    theme_classic() +
    scale_y_continuous(expand = c(0, 0), limits=c(-1,1), breaks = seq(-1, 1, by = .5)) +
    scale_fill_manual(values=c('#6E727F','#BEC1CC')) +
            theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks.length=unit(.15, "cm"),
          plot.title = element_text(hjust = 0.5, size=12),
          axis.line.x = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
          axis.text.x = element_text(colour="black", size=8, angle = 90, hjust = 1),
          axis.text.y = element_text(colour="black", size=8),
          axis.title.x = element_text(colour="black", size=8),
          axis.title.y = element_blank(),
            legend.position = "none") +
    coord_flip()


figure_out = paste0(base_dir,'/figures/PaperFig_celltype_VisDfc_cellCor_to_SSTandPVALB_AHBAGene.pdf')
figure_out
pdf(figure_out, width=3, height=3)
plot_grid(cor2pvalb, cor2sst)
dev.off()





dfc_cell2sst = ggplot(data=dfc_cor_to_sstGene, aes(x=cell, y=cor2sst)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none") +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-1, 1, by=.5),2), limits=c(-1,1))

vis_cell2sst = ggplot(data=vis_cor_to_sstGene, aes(x=cell, y=cor2sst)) +
    geom_bar(stat="identity",position="dodge") +
    coord_flip() +
            theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none") +
        scale_y_continuous(expand = c(0,.0), breaks=round(seq(-1, 1, by=.5),2), limits=c(-1,1))

figure_out = paste0(base_dir,'/figures/PaperFig_celltype_VisDfc_cellCor_to_sstGene.pdf')
pdf(figure_out, width=3, height=3)
plot_grid(dfc_cell2sst, vis_cell2sst)
dev.off()




plot_info = NULL
plot_info[['PVALB']] = list(c(0,.18,.09), c(0,.06,.03))
plot_info[['SST']]   = list(c(0,.08,.04), c(0,.08,.04))
plot_info[['In1']]   = list(c(0,.06,.03), c(0,.08,.04))
plot_info[['In3']]   = list(c(0,.12,.06), c(0,.08,.04))
plot_info[['In4']]   = list(c(0,.08,.04), c(0,.06,.03))
plot_info[['Ex1']]   = list(c(0,.24,.12), c(0,.30,.15))
plot_info[['Ex3']]   = list(c(0,.34,.17), c(0,.52,.26))
plot_info[['Ex4']]   = list(c(.08,.28,.1), c(.12,.36,.12))
plot_info[['Ex8']]   = list(c(0,.08,.04), c(.0,.10,.05))
plot_info[['Ast']]   = list(c(0,.14,.07), c(.0,.08,.04))
plot_info[['Oli']]   = list(c(0,.24,.12), c(.0,.14,.07))
plot_info[['OPC']]   = list(c(0,.04,.02), c(.0,.01,.005))
plot_info[['Mic']]   = list(c(.01,.05,.02), c(.01,.04,.015))
plot_info[['Per']]   = list(c(0,.02,.01), c(0,.02,.01))
plot_info[['End']]   = list(c(.01,.03,.01), c(0,.04,.02))

cell = 'Ex8'

figure_out = paste0(base_dir,'/figures/PaperFig_cibersort_vis_to_dfc_', cell, '.pdf')

x_lower = plot_info[[cell]][[1]][1]
x_upper = plot_info[[cell]][[1]][2]
x_break = plot_info[[cell]][[1]][3]
y_lower = plot_info[[cell]][[2]][1]
y_upper = plot_info[[cell]][[2]][2]
y_break = plot_info[[cell]][[2]][3]

plot_me = data.frame(x=as.numeric(lake_dfc_dat[lake_dfc_dat$gene == paste0('DFC_',cell), ]),
                     y=as.numeric(lake_vis_dat[lake_vis_dat$gene == paste0('VIS_',cell), ]))
summary(plot_me$x)
summary(plot_me$y)
cor.test(plot_me$x, plot_me$y)

figure_out
CairoPDF(figure_out, width=1, height=1)
print(ggplot(data=plot_me, aes(x=x, y=y)) +
        geom_point(size=.75, shape=21, stroke=.2, color='black', fill='white') +
        theme_minimal() +
        geom_smooth(method = "lm", se = FALSE, color='black', size=.5) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.ticks.length=unit(.15, "cm"),
              axis.ticks = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5, size=20),
              axis.line.y = element_line(size=0.5, linetype="solid", colour = "black"),
              axis.line.x =element_line(size=0.5, linetype="solid", colour = "black"),
              axis.text.x = element_text(colour="black",size=6),
              axis.text.y = element_text(colour="black",size=6),
              axis.title.x = element_text(colour="black",size=6),
              axis.title.y = element_text(colour="black",size=6)) +
        scale_y_continuous(expand = c(0,.0), breaks=seq(y_lower, y_upper, by=y_break), limits=c(y_lower,y_upper)) +
        scale_x_continuous(expand = c(0,.0), breaks=seq(x_lower, x_upper, by=x_break), limits=c(x_lower,x_upper)) +
        xlab(paste0('DFC ', cell, ' %')) +
        ylab(paste0('VIS ', cell, ' %')))
dev.off()








