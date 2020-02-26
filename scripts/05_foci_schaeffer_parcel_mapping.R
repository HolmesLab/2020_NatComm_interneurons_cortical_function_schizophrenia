library(tidyverse)
library(cifti)
library(gifti)
library(XML)
library(DescTools)
library(Cairo)

# function to plot parcel-wise metric files for viewing with HCP wb_view
plot_matlab = function(values, out_path, parcel_num, net_num){

    base_dir         = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'
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


# project directory
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn'

# read sample information
load(file=paste0(base_dir, '/data/ahba/donorDat_obj.Rdata'), verbose=T)
load(file=paste0(base_dir, '/data/ahba/ahba_ctx_norm.Rdata'), verbose=T)


# read sample-to-vertex projection info
sample_info = read_csv(paste0(base_dir, '/data/ahba/sample_info_vertex_mapped.csv'))
sample_dat  = sample_info[which(abs(sample_info$mm_to_surf) < 4),]


# cortical sample by normalized gene expression data frame
reg_micro_scale = as.data.frame(t(ctx_data_scale))
reg_samples = sample_dat


# stitch the left/right verticies together to match Schaeffer parcel cifti format
reg_samples$bihemi_vertex = reg_samples$vertex + 1 # cifti indices index at 0, R indexes at 1
right_ctx_idxs = intersect(grep('right', reg_samples$structure_name), which(reg_samples$top_level == 'CTX'))
reg_samples$bihemi_vertex[right_ctx_idxs] = reg_samples$bihemi_vertex[right_ctx_idxs] + 32492


# Read Schaeffer parcel info (for multiple parcel #'s, network assignments)
donor_arr =  c('9861','10021','12876','14380','15496','15697')
donor_specific_expression = NULL
for (parcels in c('200','400')){
    for (net in c('17','7')){

        # schaeffer parcellation by vertex
        schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
        schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)

        # corresponding parcel labels
        schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
        schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]

        head(schaef_cii$data)
        length(schaef_cii$data)
        sort(unique(schaef_cii$data))

        # calculate the average gene expression within each parcel
        sst_pvalb_arr = NULL
        schaeffer_mat = matrix(NA, ncol=as.numeric(parcels), nrow=ncol(reg_micro_scale))
        for (donor in donor_arr){donor_specific_expression[[donor]] = schaeffer_mat}

        # loop over each parcel, find matching samples, calculate average expression
        for (idx in 1:length(schaef_labels)){
            write(idx,'')
            # schaeffer indices
            parcel_idxs   = which(schaef_cii$data == idx)
            # foci within this parcel
            match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs)
            # data for this parcel
            match_samples = reg_samples[match_idxs,]

            # average expressino of every gene, across samples in this parcel
            schaeffer_expr = colMeans(reg_micro_scale[match_idxs,])

            # plug in values to the pre-allocated matrix
            schaeffer_mat[,idx] = schaeffer_expr

            # do the same parcel-wise averaging, but separately for each donor
            for (donor in donor_arr){
                donor_idxs     = which(as.character(reg_samples$brain) == donor)
                donor_reg_idxs = intersect(match_idxs, donor_idxs)
                expr           = colMeans(reg_micro_scale[donor_reg_idxs,])
                # average expressino of every gene, across samples in this parcel
                donor_specific_expression[[donor]][,idx] = expr
            }
        }
        # add column/row names
        schaef_out = as.data.frame(schaeffer_mat)
        rownames(schaef_out) = colnames(reg_micro_scale)
        colnames(schaef_out) = schaef_labels

        # write the averaged expression matrices for each donor
        for (donor in donor_arr){
            donor_specific_expression[[donor]] = as.data.frame(donor_specific_expression[[donor]])
            rownames(donor_specific_expression[[donor]]) = colnames(reg_micro_scale)
            colnames(donor_specific_expression[[donor]]) = schaef_labels

            donor_specific_expression[[donor]]$gene = colnames(reg_micro_scale)
            out_path = paste0(base_dir, '/data/ahba/schaeffer_donor',as.character(donor),'_p', parcels,'_',net,'Net_expr_mat.csv')
            write_csv(donor_specific_expression[[donor]], path=out_path)
        }
        sst   = as.numeric(schaef_out[rownames(schaef_out) == 'SST',])
        pvalb = as.numeric(schaef_out[rownames(schaef_out) == 'PVALB',])
        print(cor.test(sst,pvalb))

        schaef_out$gene = colnames(reg_micro_scale)
        write_csv(schaef_out, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'Net_expr_mat.csv'))
    }
}


# calculate sample coverage for each parcel
for (parcels in c('200','400')){
    net = '17'
    # schaeffer parcellation by vertex
    schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order.dscalar.nii')
    schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)

    # corresponding parcel labels
    schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
    schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]

    head(schaef_cii$data)
    length(schaef_cii$data)
    sort(unique(schaef_cii$data))

    # for each parcel, calc the # of donors with a sample present, and calc how over-represented any individual donor is
    parcel_max_brain = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    parcel_brain_compisition = matrix(NA, ncol=as.numeric(parcels), nrow=1)
    for (idx in 1:length(schaef_labels)){
        write(idx,'')
        parcel_idxs   = which(schaef_cii$data == idx) # schaeffer indices
        match_idxs    = which(reg_samples$bihemi_vertex %in% parcel_idxs) # foci within this parcel
        match_samples = reg_samples[match_idxs,] # data for this parcel

        # number of donors contributing data to this parcel
        num_brains_in_parcel = length(unique(match_samples$brain))

        # donor with maximal representation in parcel
        max_brain = names(rev(sort(table(match_samples$brain))))[1]
        max_brain_perc = length(which(match_samples$brain == max_brain))/length(match_samples$brain)

        parcel_brain_compisition[,idx] = num_brains_in_parcel
        parcel_max_brain[,idx] = max_brain_perc
    }
    colnames(parcel_brain_compisition) = schaef_labels
    colnames(parcel_max_brain) = schaef_labels
    parcel_brain_compisition = as.data.frame(parcel_brain_compisition)
    parcel_max_brain = as.data.frame(parcel_max_brain)

    write_csv(parcel_brain_compisition, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_numDonorsInEachParcel.csv'))
    write_csv(parcel_max_brain, path=paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_maxDonorPresenceInParcel.csv'))
}



# read the ROI summarized expression data from above
# plot interneuron marker distributions, summarized across scheaffer parcels
schaeffer_mat = read_csv("/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data/ahba/schaeffer400_7Net_expr_mat.csv")
pvalb_expr = as.numeric(schaeffer_mat[which(schaeffer_mat$gene == 'PVALB'),1:400])
sst_expr   = as.numeric(schaeffer_mat[which(schaeffer_mat$gene == 'SST'),1:400])
cor.test(pvalb_expr, sst_expr)


# create averaged plots
nparcel='400'
net_num='7'
plot_matlab(values=as.numeric(pvalb_expr), out_path=paste0(base_dir, '/figures/surface_plots/pvalb_expr_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')
plot_matlab(values=as.numeric(sst_expr), out_path=paste0(base_dir, '/figures/surface_plots/sst_expr_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')

donor_dist = read_csv(paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_maxDonorPresenceInParcel.csv'))
plot_matlab(values=as.numeric(donor_dist), out_path=paste0(base_dir, '/figures/surface_plots/MaxDonorPresenceInParcel_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')

donor_dist = read_csv(paste0(base_dir, '/data/ahba/schaeffer',parcels,'_',net,'_numDonorsInEachParcel.csv'))
plot_matlab(values=as.numeric(donor_dist), out_path=paste0(base_dir, '/figures/surface_plots/numDonorsInEachParcel_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')



# write the parcel-wise expression for each donor
nparcel = '400'
net_num = '7'
for (donor in donor_arr){
    path = paste0(base_dir, '/data/ahba/schaeffer_donor',as.character(donor),'_p', nparcel,'_',net_num,'Net_expr_mat.csv')
    print(path)

    donor_expr  = read_csv(path)
    donor_sst   = as.numeric(t(donor_expr[donor_expr$gene == 'SST',1:400]))
    donor_pvalb = as.numeric(t(donor_expr[donor_expr$gene == 'PVALB',1:400]))
    donor_sst_min_pvalb = donor_sst - donor_pvalb
    print(cor.test(donor_sst, donor_pvalb))

    plot_matlab(values=as.numeric(donor_sst), out_path=paste0(base_dir, '/figures/surface_plots/SST_donor',donor,'_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')
    plot_matlab(values=as.numeric(donor_pvalb), out_path=paste0(base_dir, '/figures/surface_plots/PVALB_donor',donor,'_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')
    plot_matlab(values=as.numeric(donor_sst_min_pvalb), out_path=paste0(base_dir, '/figures/surface_plots/SST_minus_PVALB_expr_donor',donor,'_schaef',nparcel,'_net',net_num,'.dscalar.nii'), parcel_num='400', net_num='7')

    all_sst_cors = cor(donor_sst, t(donor_expr[,1:400]), use='pairwise.complete')
    colnames(all_sst_cors) = donor_expr$gene
    print(length(which(all_sst_cors <= all_sst_cors[colnames(all_sst_cors) == 'PVALB']))/length(all_sst_cors))

    all_pvalb_cors = cor(donor_pvalb, t(donor_expr[,1:400]), use='pairwise.complete')
    colnames(all_pvalb_cors) = donor_expr$gene
    print(length(which(all_pvalb_cors >= all_pvalb_cors[colnames(all_pvalb_cors) == 'SST']))/length(all_pvalb_cors))
}



# create a data frame of the SST/PVALB relationship for each donor
all_cors = NULL
all_cors[['sst']] = NULL
all_cors[['pvalb']] = NULL
indiv_cor_arr = NULL
for (donor in donor_arr){
    write(donor,'')

    # donor specific cortical data
    donor_dat = reg_micro_scale[which(as.character(reg_samples$brain) == donor),]

    # get pvalb/sst expression for donor
    donor_sst   = as.numeric(t(donor_dat[,colnames(donor_dat) == 'SST']))
    donor_pvalb = as.numeric(t(donor_dat[,colnames(donor_dat) == 'PVALB']))


    print(cor.test(donor_sst, donor_pvalb))
    print(cor.test(donor_sst, donor_pvalb, method='spearman'))

    # pearson correlation for each donor
    all_sst_cors        = cor(donor_sst, donor_dat, use='pairwise.complete')
    all_cors[['sst']]   = rbind(all_cors[['sst']], all_sst_cors)
    all_pvalb_cors      = cor(donor_pvalb, donor_dat, use='pairwise.complete')
    all_cors[['pvalb']] = rbind(all_cors[['pvalb']], all_pvalb_cors)

    all_donor_cors = cor(donor_dat)
    cor_arr   = all_donor_cors[upper.tri(all_donor_cors)]
    perm_pval = length(which(cor_arr <= cor(donor_sst, donor_pvalb)))/length(cor_arr)

    sp_all_donor_cors = cor(donor_dat, method='spearman')
    sp_cor_arr   = sp_all_donor_cors[upper.tri(sp_all_donor_cors)]
    spear_perm_pval = length(which(sp_cor_arr <= cor(donor_sst, donor_pvalb, method='spearman')))/length(sp_cor_arr)

    sub_sst_pvalb = cor.test(donor_sst, donor_pvalb)
    out_row       = data.frame(type='pearson', perm_pval=perm_pval, spear_perm_pval=spear_perm_pval, cor=sub_sst_pvalb$estimate, donor=donor, ci95low = sub_sst_pvalb$conf.int[1], ci95high = sub_sst_pvalb$conf.int[2])
    indiv_cor_arr = rbind(indiv_cor_arr, out_row)

    sub_sst_pvalb = as.data.frame(t(SpearmanRho(donor_sst, donor_pvalb, conf.level = 0.95)))
    out_row       = data.frame(type='spearman', perm_pval=perm_pval, spear_perm_pval=spear_perm_pval, cor=sub_sst_pvalb$rho, donor=donor, ci95low = sub_sst_pvalb$lwr.ci, ci95high = sub_sst_pvalb$upr.ci)
    indiv_cor_arr = rbind(indiv_cor_arr, out_row)
}
collapsed_sst_cors   = data.frame(cors = colMeans(all_cors[['sst']]))
collapsed_pvalb_cors = data.frame(cors = colMeans(all_cors[['pvalb']]))



out_file = paste0(base_dir, '/figures/PaperFig_individual_CortCors_sst_pvalb.pdf')
out_file
CairoPDF(out_file, width=2, height=4)
indiv_cor_arr$color='white'
indiv_cor_arr$idx = rev(c(1,2,4,5,7,8,10,11,13,14,16,17)) #1:length(indiv_cor_arr$idx)
ggplot(indiv_cor_arr, aes(x=idx, y=cor)) +
    geom_point(size=3) +
    coord_flip() +
    theme_minimal() +
    expand_limits(y = c(-0.75, 0)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(-.6,0,by=0.2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-.75,0,by=0.25)) +
    geom_linerange(aes(ymin=indiv_cor_arr$ci95low, ymax=indiv_cor_arr$ci95high), position=position_dodge(0.4), size=.5)
dev.off()



# Control analysis that isn't in the actual paper
# provides a good sanity check of the SST/PVALB correlation effects for interested parties

# basically, we split groups in two and show that between-groups parcel-wise correlation of SST/PVALB 
# rules out any potential within-subjects effects as a driver of the correlation.  
n_donor = 6
x  = 1:n_donor
x1 = combn(x, n_donor/2) #how many ways can we take half the elements to form the 1st group
NC = NCOL(x1)
x2 = x1[, NC:1] # simplified way to generate the complementary groups that include values not in x1
grp1 = t(x1[,1:(NC/2)]) # We only need half of the rows, the 2nd half containing the same set in reverse order
grp2 = t(x2[,1:(NC/2)])
all_comb = cbind(grp1, grp2)

parcels = '400'
net = '7'
# schaeffer parcellation by vertex
schaeffer  = paste0(base_dir, '/data/Schaefer/Schaefer2018_', parcels,'Parcels_', net,'Networks_order.dscalar.nii')
schaef_cii = read_cifti(schaeffer, drop_data = TRUE, trans_data = TRUE)

# corresponding parcel labels
schaeffer_labels = read_csv(col_names=F, paste0(base_dir, '/data/Schaefer/Schaefer2018_',parcels,'Parcels_',net,'Networks_order_info.txt'))
schaef_labels    = schaeffer_labels$X1[grep('Network', schaeffer_labels$X1)]

# for every group-wise combination
group_expr_list = NULL
for (r in 1:nrow(all_comb)){
    print(r)
    row = all_comb[r,]
    grp1 = row[1:3]
    grp2 = row[4:6]

    # sample idxs from donors in the 1st/2nd group
    grp1_idxs = which(as.character(reg_samples$brain) %in% donor_arr[grp1])
    grp2_idxs = which(as.character(reg_samples$brain) %in% donor_arr[grp2])

    # pre_allocate output
    grp_expr = NULL
    grp_expr[['grp1']] = matrix(NA, ncol=as.numeric(parcels), nrow=ncol(reg_micro_scale))
    grp_expr[['grp2']] = matrix(NA, ncol=as.numeric(parcels), nrow=ncol(reg_micro_scale))

    # For each group. Calculate the average gene expression within each parcel
    for (idx in 1:(length(schaef_labels)) ){
        # identify samples for each parcel/donor/
        # schaeffer indices
        parcel_idxs = which(schaef_cii$data == idx)
        # foci within this parcel
        match_idxs  = which(reg_samples$bihemi_vertex %in% parcel_idxs)
        grp1_use    = intersect(grp1_idxs, match_idxs)
        grp2_use    = intersect(grp2_idxs, match_idxs)

        # average expressino of every gene, across samples in this parcel
        schaeffer1_expr = colMeans(reg_micro_scale[grp1_use,])
        # average expressino of every gene, across samples in this parcel
        schaeffer2_expr = colMeans(reg_micro_scale[grp2_use,])

        # plug in values to the pre-allocated matrix
        grp_expr[['grp1']][,idx] = schaeffer1_expr
        grp_expr[['grp2']][,idx] = schaeffer2_expr
    }
    rownames(grp_expr[['grp1']]) = colnames(reg_micro_scale)
    rownames(grp_expr[['grp2']]) = colnames(reg_micro_scale)

    # compare group_1 SST to group_2 PVALB
    grp1_sst   = grp_expr[['grp1']][rownames(grp_expr[['grp1']]) == 'SST',1:200]
    grp2_pvalb = grp_expr[['grp2']][rownames(grp_expr[['grp2']]) == 'PVALB',1:200]
    corA = cor.test(grp1_sst, grp2_pvalb)
    all_corsA = cor(grp1_sst, t(grp_expr[['grp2']][,1:200]), use='pairwise.complete')
    length(which(all_corsA  < corA$estimate))

    # compare group_2 SST to group_1 PVALB
    grp2_sst = grp_expr[['grp2']][rownames(grp_expr[['grp2']]) == 'SST',1:200]
    grp1_pvalb = grp_expr[['grp1']][rownames(grp_expr[['grp1']]) == 'PVALB',1:200]
    corB = cor.test(grp2_sst, grp1_pvalb)
    all_corsB = cor(grp2_sst, t(grp_expr[['grp1']][,1:200]), use='pairwise.complete')
    length(which(all_corsB  < corB$estimate))

    group_expr_list[[r]] = grp_expr
}
# "group_expr_list" contains all the averaged parcel-wise expression for each group split



sst_pvalb_cors  = NULL
sst_global_cors = NULL
pvalb_global_cors = NULL
for (r in 1:nrow(all_comb)){

    print(r)
    # compare
    cur_exprA  = group_expr_list[[r]][['grp1']]
    cur_exprB  = group_expr_list[[r]][['grp2']]

    grp1_sst   = cur_exprA[rownames(cur_exprA) == 'SST',1:200]
    grp1_pvalb = cur_exprA[rownames(cur_exprA) == 'PVALB',1:200]

    grp2_sst   = cur_exprB[rownames(cur_exprB) == 'SST',1:200]
    grp2_pvalb = cur_exprB[rownames(cur_exprB) == 'PVALB',1:200]

    # SST/PVALB pairwise correlation between groups
    corA = cor.test(grp1_sst, grp2_pvalb)
    print(corA$estimate)
    corB = cor.test(grp2_sst, grp1_pvalb)
    print(corB$estimate)

    sst_pvalb_cors = c(sst_pvalb_cors, corA$estimate, corB$estimate)

    # left hemi sst for group 1/2
    all_grp1_sst = cor(grp1_sst, t(cur_exprB[,1:200]), use='pairwise.complete')
    all_grp2_sst = cor(grp2_sst, t(cur_exprA[,1:200]), use='pairwise.complete')

    # left hemi pvalb for group 1/2
    all_grp1_pvalb = cor(grp1_pvalb, t(cur_exprB[,1:200]), use='pairwise.complete')
    all_grp2_pvalb = cor(grp2_pvalb, t(cur_exprA[,1:200]), use='pairwise.complete')

    sst_global_cors = rbind(sst_global_cors, all_grp1_sst)
    sst_global_cors = rbind(sst_global_cors, all_grp2_sst)

    pvalb_global_cors = rbind(pvalb_global_cors, all_grp1_pvalb)
    pvalb_global_cors = rbind(pvalb_global_cors, all_grp2_pvalb)

    print(length(which(all_grp1_sst <= corA$estimate)))
    print(length(which(all_grp2_sst <= corB$estimate)))
}

avg_pvalb_cors_df = data.frame(avg_cor=colMeans(pvalb_global_cors), gene=names(colMeans(pvalb_global_cors)))
avg_pvalb_cors_df = avg_pvalb_cors_df[rownames(avg_pvalb_cors_df) != 'PVALB',]
avg_pvalb_cors_df$fill = 'null'
avg_pvalb_cors_df$fill[avg_pvalb_cors_df$gene == 'SST'] = 'SST'
avg_pvalb_cors_df$size = 1
avg_pvalb_cors_df$size[rownames(avg_pvalb_cors_df) == 'SST'] = 1.5
avg_pvalb_cors_df = avg_pvalb_cors_df[order(avg_pvalb_cors_df$avg_cor),]
avg_pvalb_cors_df$idx = 1:nrow(avg_pvalb_cors_df)


# benchmark SST against all other gene-wise correlations to PVALB
length(which(avg_pvalb_cors_df$avg_cor <= avg_pvalb_cors_df$avg_cor[avg_pvalb_cors_df$gene == 'SST']))

# position of the SST/PVALB given the distribution
pop_mean = mean(avg_pvalb_cors_df$avg_cor[avg_pvalb_cors_df$gene != 'SST'])
pop_sd   = sd(avg_pvalb_cors_df$avg_cor[avg_pvalb_cors_df$gene != 'SST'])
sst_pvalb_cor = avg_pvalb_cors_df$avg_cor[avg_pvalb_cors_df$gene == 'SST']


CairoPDF('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_split_group_SstPvalb_PVALB_reference.pdf', family='Arial', width=5, height=3)
sst_value = avg_pvalb_cors_df$avg_cor[which(avg_pvalb_cors_df$gene == 'SST')]
pop_mean  = mean(avg_pvalb_cors_df$avg_cor)
pop_sd  = mean(avg_pvalb_cors_df$avg_cor)

ggplot(avg_pvalb_cors_df, aes(avg_cor)) +
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
    geom_vline(xintercept = sst_value, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()


avg_sst_cors_df = data.frame(avg_cor=colMeans(sst_global_cors), gene=names(colMeans(sst_global_cors)))
avg_sst_cors_df = avg_sst_cors_df[rownames(avg_sst_cors_df) != 'SST',]
avg_sst_cors_df = avg_sst_cors_df[order(avg_sst_cors_df$avg_cor),]
avg_sst_cors_df$fill = 'null'
avg_sst_cors_df$fill[rownames(avg_sst_cors_df) == 'PVALB'] = 'PVALB'
avg_sst_cors_df$size = 1
avg_sst_cors_df$size[rownames(avg_sst_cors_df) == 'PVALB'] = 1.5
avg_sst_cors_df = avg_sst_cors_df[rownames(avg_sst_cors_df) != 'SST',]
avg_sst_cors_df = avg_sst_cors_df[order(avg_sst_cors_df$avg_cor),]
avg_sst_cors_df$idx = 1:nrow(avg_sst_cors_df)

# benchmark PVALB against all other gene-wise correlations to SST
length(which(avg_sst_cors_df$avg_cor <= avg_sst_cors_df$avg_cor[avg_sst_cors_df$gene == 'PVALB']))

# position of SST/PVALB given the distribution
pop_mean = mean(avg_sst_cors_df$avg_cor[avg_pvalb_cors_df$gene != 'PVALB'])
pop_sd   = sd(avg_sst_cors_df$avg_cor[avg_sst_cors_df$gene != 'PVALB'])
sst_pvalb_cor = avg_sst_cors_df$avg_cor[avg_sst_cors_df$gene == 'PVALB']


CairoPDF('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/figures/PaperFig_split_group_SstPvalb_SST_reference.pdf', family='Arial', width=5, height=3)
pvalb_value = avg_sst_cors_df$avg_cor[which(avg_sst_cors_df$gene == 'PVALB')]
pop_mean  = mean(avg_sst_cors_df$avg_cor)
pop_sd  = sd(avg_sst_cors_df$avg_cor)

ggplot(avg_sst_cors_df, aes(avg_cor)) +
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
    geom_vline(xintercept = pvalb_value, color='black', size=1, linetype='longdash') +
    geom_vline(xintercept = pop_mean, color='black', size=1, linetype='longdash')
dev.off()




out_file = paste0(base_dir, '/figures/group_split_pvalb_dist.pdf')
out_file
CairoPDF(out_file, width=4, height=2)
ggplot(avg_pvalb_cors_df, aes(x=idx, y=avg_cor, color=fill, size=size)) +
    geom_point(stroke=0) +
    theme_minimal() +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    expand_limits(y=c(-1,1)) +
    #scale_x_continuous(expand = c(0.25,.25) ) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-1, 1, by = .5))
dev.off()


avg_sst_cors_df = data.frame(avg_cor=colMeans(sst_global_cors), gene=names(colMeans(sst_global_cors)))
avg_sst_cors_df = avg_sst_cors_df[order(avg_sst_cors_df$avg_cor),]
avg_sst_cors_df$fill = 'null'
avg_sst_cors_df$fill[rownames(avg_sst_cors_df) == 'PVALB'] = 'PVALB'
avg_sst_cors_df$size = 1
avg_sst_cors_df$size[rownames(avg_sst_cors_df) == 'PVALB'] = 1.5
avg_sst_cors_df = avg_sst_cors_df[rownames(avg_sst_cors_df) != 'SST',]
avg_sst_cors_df = avg_sst_cors_df[order(avg_sst_cors_df$avg_cor),]
avg_sst_cors_df$idx = 1:nrow(avg_sst_cors_df)

out_file = paste0(base_dir, '/figures/group_split_sst_dist.pdf')
out_file
CairoPDF(out_file, width=4, height=2)
ggplot(avg_sst_cors_df, aes(x=idx, y=avg_cor, color=fill, size=size)) +
    geom_point(stroke=0) +
    theme_minimal() +
    scale_color_manual(values=c("#999999", "#E69F00")) +
    expand_limits(y=c(-1,1)) +
    #scale_x_continuous(expand = c(0.25,.25) ) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(-1, 1, by = .5))
dev.off()


avg_sst_cors = data.frame(avg_cor=colMeans(sst_global_cors))
avg_sst_cors$fill = 'null'
avg_sst_cors$fill[rownames(avg_sst_cors) == 'PVALB'] = 'PVALB'
avg_sst_cors = avg_sst_cors[order(avg_sst_cors$avg_cor),]
avg_sst_cors$idx = 1:nrow(avg_sst_cors)



# Foci mapping
# read previously created file with vertex and sst/pvalb information

# source project utility library
source(paste0(base_dir, '/scripts/R_function_library.R'))

# read vertex mapping file
vert.dir = paste0(base_dir, '/data/ahba')
vert_ctx = sample_dat


# Mean-normalized cortical tissue samples
reg = 'CTX'
ctx_micro_scale = reg_micro_scale
rownames(ctx_micro_scale) = vert_ctx$well_id

# get difference betw SST-PVALB, throw it in the cortical df
sst_v_pvalb    = ctx_micro_scale$SST - ctx_micro_scale$PVALB
vert_ctx$sst_v_pvalb = sst_v_pvalb
vert_ctx$sst   = ctx_micro_scale$SST
vert_ctx$pvalb = ctx_micro_scale$PVALB


# donor names
sub_list = c('9861', '10021', '12876', '14380', '15496', '15697')

wb_path = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'


# folder to save surface plots into
# lh9861-164
hgnc_A = 'SST' # sst
hgnc_B = 'PVALB' # pvalb
name_base = 'SST_minus_PVALB'
out_dir = paste0(base_dir, '/figures/surface_plots/')
for (hemi in c('lh','rh')){
    for ( donor in sub_list ){
        if (hemi=='lh'){
            name_string   = paste0(name_base,'_LH')
            vert_ctx_hemi = vert_ctx[intersect(which(vert_ctx$brain == as.numeric(donor)), grep('left', vert_ctx$structure_name)),]
            ctx_expr      = ctx_micro_scale[which(rownames(ctx_micro_scale) %in% vert_ctx_hemi$well_id),]
            ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'
            projection_structure = 'CORTEX_LEFT'

        } else if (hemi=='rh') {
            name_string   = paste0(name_base,'_RH')
            vert_ctx_hemi = vert_ctx[intersect(which(vert_ctx$brain == as.numeric(donor)), grep('right', vert_ctx$structure_name)),]
            ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'
            ctx_expr      = ctx_micro_scale[which(rownames(ctx_micro_scale) %in% vert_ctx_hemi$well_id),]
            projection_structure = 'CORTEX_RIGHT'

            if (!donor %in% c('9861','10021')){
                next
            }
        }
        # foci file name
        out_name = paste(out_dir, 'PaperFoci_donor', donor, '_', name_string, '_', projection_structure, '.foci', sep='')
        print(out_name)

        # vertices and MNI coords
        mni_coords = vert_ctx_hemi[c('mni_x','mni_y','mni_z','R','A','S')]
        vertices   = vert_ctx_hemi$vertex

        # convert expression values to rgb
        vals_A = ctx_expr[,which(colnames(ctx_expr) == hgnc_A)]
        vals_B = ctx_expr[,which(colnames(ctx_expr) == hgnc_B)]
        plot_vals = vals_A - vals_B

        mypal      = colorRampPalette( c( '#4C5EAA','#7D94CB','#B8CCEA', "white", '#EED0BF','#E9795F','#BD232D') )( 100 )
        rgb_matrix = t(col2rgb(map2color(plot_vals, mypal, limits=c(-3,3))))/256
        rgba_color = cbind(rgb_matrix, rep(1,dim(rgb_matrix)[1]))
        rgba_color = as.data.frame(rgba_color)
        colnames(rgba_color) = c('r', 'g', 'b', 'a')
        rgba_color = as.data.frame(rgba_color)

        # do the plotting
        plot_foci(mni.coords=mni_coords,
                  vertices=vertices,
                  foci.names=vert_ctx_hemi$well_id,
                  ref.surface=ref_surface,
                  rgba.color=rgba_color,
                  wb.path=wb_path,
                  projection.structure=projection_structure,
                  out.name=out_name)
  }
}
# color scale for figure 1
z=matrix(1:100,nrow=1)
x=1
y=seq(-3,3,len=100)
image(x,y,z,col=mypal,axes=FALSE,xlab="",ylab="")




# create individual expression plots for each gene of interest
for (gene in c('sst','pvalb')){
  for (hemi in c('lh','rh')){
    for ( donor in sub_list ){
      if (hemi=='lh'){
        name_string   = paste0(gene, '_LH')
        vert_ctx_hemi = vert_ctx[intersect(which(vert_ctx$brain == as.numeric(donor)), grep('left', vert_ctx$structure_name)),]
        ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.L.inflated_MSMAll.32k_fs_LR.surf.gii'
        projection_structure = 'CORTEX_LEFT'

      } else if (hemi=='rh') {
        name_string   =  paste0(gene, '_RH')
        vert_ctx_hemi = vert_ctx[intersect(which(vert_ctx$brain == as.numeric(donor)), grep('right', vert_ctx$structure_name)),]
        ref_surface   = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/HCP_S1200_GroupAvg_v1/S1200.R.inflated_MSMAll.32k_fs_LR.surf.gii'
        projection_structure = 'CORTEX_RIGHT'
        if (!donor %in% c('9861','10021')){
          next
        }
      }
      if (gene == 'sst'){
        mypal = colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D') )( 100 )
      } else if (gene == 'pvalb'){
        mypal = colorRampPalette( c("white", '#EBF0B5', '#46B6C2', '#283B8E')  )( 100 )
      }
      out_name = paste(out_dir, 'PaperFoci_donor', donor, '_', name_string, '_', projection_structure, '.foci', sep='')
      print(out_name)
      # subset LH and donor sample information

      # vertices and MNI coords
      mni_coords = vert_ctx_hemi[c('mni_x','mni_y','mni_z','R','A','S')]
      vertices   = vert_ctx_hemi$vertex

      # convert expression values to rgb
      rgb_matrix = t(col2rgb(map2color(vert_ctx_hemi[[gene]], mypal, limits=c(-2.0, 2.0))))/256
      rgba_color = cbind(rgb_matrix, rep(1,dim(rgb_matrix)[1]))
      rgba_color = as.data.frame(rgba_color)
      colnames(rgba_color) = c('r', 'g', 'b', 'a')
      rgba_color = as.data.frame(rgba_color)

      # Do the plotting
      plot_foci(mni.coords=mni_coords,
                  vertices=vertices,
                  foci.names=vert_ctx_hemi$well_id,
                  ref.surface=ref_surface,
                  rgba.color=rgba_color,
                  wb.path=wb_path,
                  projection.structure=projection_structure,
                  out.name=out_name)
    }
  }
}

# create color scale for each gene
for (gene in c('sst','pvalb')){
  if (gene == 'sst'){
    mypal <- colorRampPalette( c("white", '#FCEDA9', '#EF924F', '#AE232D') )( 100 )
  } else if (gene == 'pvalb'){
    mypal <- colorRampPalette( c("white", '#EBF0B5', '#46B6C2', '#283B8E') )( 100 )
  }
  # color scale for figure 1
  pdf(paste0(base_dir, '/figures/', gene, '_supp_fig.pdf'))
  z=matrix(1:100,nrow=1)
  x=1
  y=seq(-3,3,len=100)
  print(image(x,y,z,col=mypal,axes=FALSE,xlab="",ylab=""))
  dev.off()
}








