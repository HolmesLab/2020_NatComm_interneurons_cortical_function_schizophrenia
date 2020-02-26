library(Cairo)
library(tidyverse)


# Modify these filepaths for your local directory structure

# enter the abin directory for local install of AFNI
afni_dir = '/gpfs/milgram/apps/hpc.rhel7/software/AFNI/2018.08.28' 

# enter the directory containing the script repository
base_dir = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn' 


# read previously processed AHBA data
load_me = paste0(base_dir, '/data/ahba/donorDat_obj.Rdata')
load_me
load(file=load_me, verbose=T)


# Source function library for this project
source(paste(base_dir, '/scripts/R_function_library.R', sep = ''))


# Identify which samples overlap with Yeo/Choi Atlases, using AHBA provided MNI locations
atlas_dir  = paste0(base_dir, '/data/atlas/')
donor_nums = unique(donorDat$samples$brain)


# z-transform striatal samples for each donor
str_samples = NULL
str_micros_scale = NULL
for (donor in donor_nums ){
    write(donor,'')
    str_donor_idxs = intersect(which(donorDat$samples$top_level == 'STR'), which(donorDat$samples$brain == donor))

    str_samp    = donorDat$samples[str_donor_idxs,]
    str_samples = rbind(str_samples, str_samp)
    str_micro   = donorDat$micro[,str_donor_idxs]

    str_scale   = apply(str_micro, 1, scale)
    str_micros_scale = cbind(str_micros_scale, as.matrix(t(str_scale)))
}
str_micros_scale = as.data.frame(str_micros_scale)
str_samples = str_samples[unique(colnames(str_samples))]

# SST and PVALB striatal expression
probes = donor_dat$probes_unique
sst_expr    = as.numeric(str_micros_scale[rownames(str_micros_scale) == 'SST',])
pvalb_expr  = as.numeric(str_micros_scale[rownames(str_micros_scale) == 'PVALB',])

sst_v_pvalb = sst_expr - pvalb_expr
cor.test(pvalb_expr, sst_expr)


# Chris Gorgolewski MNI mapping ( https://github.com/chrisgorgo/alleninf/tree/master/alleninf/data )
mni_corrected = read_csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/alleninf/alleninf/data/corrected_mni_coordinates.csv')
str_samples$idx = 1:nrow(str_samples)
str_samples   = merge(x=str_samples, y=mni_corrected, 'well_id')
str_samples   = str_samples[order(str_samples$idx),]

# call function to see if MNI coordinates of each sample overlaps with Choi 2012 atlas
atlas   = paste0(atlas_dir, 'Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz')
striat7_orig = striat_query(atlas=atlas,
                        MNI_coords=str_samples,
                        mni_x = 'mni_x',
                        mni_y = 'mni_y',
                        mni_z = 'mni_z',
                        afni.dir = afni_dir)

# Gorgolewski coordinates
atlas   = paste0(atlas_dir, 'Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz')
striat7_corrected = striat_query(atlas=atlas,
                        MNI_coords=str_samples,
                        mni_x = 'corrected_mni_x',
                        mni_y = 'corrected_mni_y',
                        mni_z = 'corrected_mni_z',
                        afni.dir = afni_dir)

# add columns for functional network assingment/gene expression
choi_names = c('none','visual','sommot','dorsattn','ventattn','limbic','control','default')
table(choi_names[striat7_orig+1], choi_names[striat7_corrected+1])


# Results were consistent across AHBA and Gorgolewski's corrected i,j,k coordinates. 
# We use the AHBA provided values
str_samples$choi_net = striat7_orig
#str_samples$choi_net = striat7_corrected
str_samples$sst_v_pvalb = as.numeric(sst_v_pvalb)
str_samples$sst      = as.numeric(sst_expr)
str_samples$pvalb    = as.numeric(pvalb_expr)


# test for differential SST/PVALB expression in limbic striatal
# limbic t-test
t.test(str_samples$sst[which(str_samples$choi_net == 5)], str_samples$pvalb[which(str_samples$choi_net == 5)], paired=T)
wilcox.test(str_samples$sst[which(str_samples$choi_net == 5)], str_samples$pvalb[which(str_samples$choi_net == 5)], paired=T)
str_samples %>% filter(choi_net==5) %>% summarize(SST_mean=mean(sst),SST_sd=sd(sst),PVALB_mean=mean(pvalb),PVALB_sd=sd(pvalb))


# test for differential SST/PVALB expression in som/mot striatal
# sommot t-test
t.test(str_samples$sst[which(str_samples$choi_net == 2)], str_samples$pvalb[which(str_samples$choi_net == 2)], paired=T)
wilcox.test(str_samples$sst[which(str_samples$choi_net == 2)], str_samples$pvalb[which(str_samples$choi_net == 2)], paired=T)
str_samples %>% filter(choi_net==2) %>% summarize(SST_mean=mean(sst),SST_sd=sd(sst),PVALB_mean=mean(pvalb),PVALB_sd=sd(pvalb))


choi_names = c('none','visual','sommot','dorsattn','ventattn','limbic','control','default')
str_samples$net_names = choi_names[str_samples$choi_net+1]

# reorganize for plotting
str_int = str_samples %>%
            group_by(net_names) %>%
            summarize(Median=median(sst),Min=min(sst),Max=max(sst),Q1=quantile(sst,.25),Q3=quantile(sst,.75),gene='sst',SSTmPVALB=median(sst)-median(pvalb))
pvalb_int = str_samples %>%
            group_by(net_names) %>%
            summarize(Median=median(pvalb),Min=min(pvalb),Max=max(pvalb),Q1=quantile(pvalb,.25),Q3=quantile(pvalb,.75),gene='pvalb',SSTmPVALB=median(sst)-median(pvalb))

plot_data = rbind(str_int, pvalb_int)
plot_data = plot_data[plot_data$net_names %in% choi_names[2:8],]
plot_data$region = factor(plot_data$net_names, levels=c('limbic','default','control','ventattn','sommot'))
plot_data$gene   = factor(plot_data$gene, levels=c('sst', 'pvalb'))


# Plot SST/PVALB expression across striatal subregions
CairoPDF(paste0(base_dir,'/figures/PaperFigure_choi7_striatum_sst_pvalb.pdf'), height=2.5, width=4)
ggplot(data=plot_data, aes(x=region, y=Median, color=gene)) +
       geom_pointrange(aes(ymin=plot_data$Q1, ymax=plot_data$Q3), position=position_dodge(0.4), size=.8) +
       geom_linerange(aes(ymin=plot_data$Min, ymax=plot_data$Max), position=position_dodge(0.4), size=.2) +
       theme_minimal() +
        scale_color_manual(values=c('#F15854','#5DA5DA')) +
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
       scale_y_continuous(expand = c(0,0), limits=c(-4,4.03)) +
       xlab("") +
       ylab("Relative SST-PVALB expression") +
       ggtitle('Choi Atlas Striatum')
dev.off()



# 9-network - Thalamus
thal_atlas = paste0(atlas_dir, 'MGH_WTA_July2016_masked_1mm.nii')


# mean-normalize thalamic atlas' for each donor
thal_samples = NULL
thal_micros_scale = NULL
for (donor in donor_nums ){
    write(donor,'')
    thal_idxs   = which(donorDat$samples$top_level == 'THAL')
    dorsal_thal = which(!donorDat$samples$structure_acronym %in% c('R','ZI'))

    thal_idxs   = intersect(thal_idxs, dorsal_thal)
    donor_idxs  = which(donorDat$samples$brain == donor)

    thal_samp    = donorDat$samples[intersect(thal_idxs, donor_idxs),]
    thal_samples = rbind(thal_samples, thal_samp)
    thal_micro   = donorDat$micro[,intersect(thal_idxs, donor_idxs)]

    thal_scale   = apply(thal_micro, 1, scale)
    thal_micros_scale = cbind(thal_micros_scale, as.matrix(t(thal_scale)))
}
thal_micros_scale = as.data.frame(thal_micros_scale)

probes      = donor_dat$probes_unique
sst_expr    = as.numeric(thal_micros_scale[rownames(thal_micros_scale) == 'SST',])
pvalb_expr  = as.numeric(thal_micros_scale[rownames(thal_micros_scale) == 'PVALB',])
sst_v_pvalb = sst_expr - pvalb_expr
cor.test(pvalb_expr, sst_expr)


# Chris Gorgo coordinates
mni_corrected = read_csv('/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/external/alleninf/alleninf/data/corrected_mni_coordinates.csv')
thal_samples$idx = 1:nrow(thal_samples)
thal_sample_new = merge(x=thal_samples, y=mni_corrected, 'well_id')
thal_sample_new = thal_sample_new[order(thal_sample_new$idx),]
thal_sample_new = thal_sample_new[unique(colnames(thal_sample_new))]


# Gorgolewski MNI
thal_corrected = striat_query(atlas=thal_atlas,
                        MNI_coords=thal_sample_new,
                        mni_x = 'corrected_mni_x',
                        mni_y = 'corrected_mni_y',
                        mni_z = 'corrected_mni_z',
                        afni.dir = afni_dir)

# AHBA MNI
thal_orig = striat_query(atlas=thal_atlas,
                    MNI_coords=thal_sample_new,
                    mni_x = 'mni_x',
                    mni_y = 'mni_y',
                    mni_z = 'mni_z',
                    afni.dir = afni_dir)


# query the thalamic atlas
#thal_sample_new$thal_net = thal_corrected
thal_sample_new$thal_net = thal_orig
thal_sample_new$sst_v_pvalb = as.numeric(sst_v_pvalb)
thal_sample_new$sst      = as.numeric(sst_expr)
thal_sample_new$pvalb    = as.numeric(pvalb_expr)


# convert ints to net names
thal_names = c('none','default','cingulopercular','sommot','frontoparietal','latOcc','medOcc','medTemp','Temp','dorsAttn')
thal_sample_new$thal_net = thal_names[thal_sample_new$thal_net+1]
use_regions = table(thal_sample_new$thal_net)[table(thal_sample_new$thal_net) > 2]


table(thal_corrected, thal_sample_new$structure_acronym)



t.test(thal_sample_new$sst[thal_sample_new$structure_acronym == 'DTLv'], thal_sample_new$pvalb[thal_sample_new$structure_acronym == 'DTLv'], paired=T)
wilcox.test(thal_sample_new$sst[thal_sample_new$structure_acronym == 'DTLv'], thal_sample_new$pvalb[thal_sample_new$structure_acronym == 'DTLv'], paired=T)


# default t-test
t.test(thal_sample_new$sst[thal_sample_new$thal_net == 'default'], thal_sample_new$pvalb[thal_sample_new$thal_net == 'default'], paired=T)
wilcox.test(thal_sample_new$sst[thal_sample_new$thal_net == 'default'],thal_sample_new$pvalb[thal_sample_new$thal_net == 'default'],paired=TRUE)

# sommot t-test
t.test(thal_sample_new$sst[thal_sample_new$thal_net == 'sommot'], thal_sample_new$pvalb[thal_sample_new$thal_net == 'sommot'], paired=T)
wilcox.test(thal_sample_new$sst[thal_sample_new$thal_net == 'sommot'],thal_sample_new$pvalb[thal_sample_new$thal_net == 'sommot'],paired=TRUE)



# get data ready for plotting
sst_thal_int = thal_sample_new %>%
            group_by(thal_net) %>%
            summarize(Median=median(sst),Min=min(sst),Max=max(sst),Q1=quantile(sst,.25),Q3=quantile(sst,.75),gene='sst',SSTmPVALB=median(sst)-median(pvalb))
pvalb_thal_int = thal_sample_new %>%
            group_by(thal_net) %>%
            summarize(Median=median(pvalb),Min=min(pvalb),Max=max(pvalb),Q1=quantile(pvalb,.25),Q3=quantile(pvalb,.75),gene='pvalb',SSTmPVALB=median(sst)-median(pvalb))

plot_data = rbind(sst_thal_int, pvalb_thal_int)
use_regions = names(use_regions)[names(use_regions) != 'none']
plot_data = plot_data[plot_data$thal_net %in% use_regions,]
plot_data$gene      = factor(plot_data$gene, levels=c('sst', 'pvalb'))
plot_data$region    = factor(plot_data$thal_net, levels=c('default','cingulopercular','frontoparietal','sommot','dorsAttn','latOcc','medOcc'))



CairoPDF(paste0(base_dir,'/figures/PaperFig_thalamus_sst_pvalb.pdf'), height=2.5, width=3.5)
ggplot(data=plot_data, aes(x=region, y=Median, color=gene)) +
  geom_pointrange(aes(ymin=plot_data$Q1, ymax=plot_data$Q3), position=position_dodge(0.4), size=.8) +
  geom_linerange(aes(ymin=plot_data$Min, ymax=plot_data$Max), position=position_dodge(0.4), size=.2) +
  theme_minimal() + scale_color_manual(values=c('#F15854','#5DA5DA')) +
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
  ylab("Relative SST-PVALB expression") +
  ggtitle('Thalamus')
dev.off()



