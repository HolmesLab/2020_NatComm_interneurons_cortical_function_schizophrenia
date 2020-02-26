# 2020_NatComm_interneurons_and_cortical_function
Code repository for "Transcriptional and imaging-genetic association of cortical interneurons, brain function, and schizophrenia risk"


**01_preprocess_ahba.R**.  
    - Reads AHBA data, gets rid of:  
        (a) probes with no entrez id.  
        (b) probes that are not expressed on average in at least 20% of cortical samples.  
    - Assigns top-level categories to AHBA ontology annotations (e.g. mediodorsal thalamus >> THAL).  
    - Saves a csv with sample information from each donor.  
    - If multiple probes are present, keep the one with the MaxMean (CollapseRows function).  
    - Saves 'donorDat_obj' with raw/subsetted AHBA data for each donor.  


**02_project_freesurfer.bash**.  
    - Raw Freesurfer AHBA surfaces downloaded from Romero and colleagues (https://www.repository.cam.ac.uk/handle/1810/265272).   
    - Calculates midthickness between pial/white.  
    - Uses HCP templates to transform individual mesh to to 32k FSLR template.  
    - Applies linear transform to move surface shape files from freesurfer back to individual coordinate reference space.  
    - Final output is a midthickness surface, deformed to native subject geometry, with aligned 32kFSLR mesh vertices.  


**03_map_ahba_ctx_to_surface.py**
    - Use the IJK coordinates of each AHBA donor tissue sample to map to the closest midthickness vertex.  


**04_ahba_sst_v_pvalb.R**
    - Main script for characterizing cortical SST v PVALB correlations in cortex.  
    - z-transforms cortical expression data separately for each donor.  
    - sample-wise correlation and scatterplot of SST and PVALB across cortex (Figure 1e).  
    - write csv files of the cortical spatial correlation of every gene to sst and pvalb.  
    - correlation of SST and PVALB in each subcortical region examined.  
    - produce SST-PVALB correlation plots for each region.  
    - compare subcortical SST+PVALB to ground truth rodent cell densities.  


**04b_cortico_striato_thalamic.R**
    - map volumetric locations of striatal/thalamic atlas of Choi et al. (2012) and Hwang et al. (2017). Both with AHBA and Chris Gorgo MNI coordinates.  
    - plot SST/PVALB for each functional defined subregion.  
    - t.test and wilcoxin test for significant differences in normalized SST and PVALB for limbic/sommot regions.  


**05_foci_schaeffer_parcel_mapping.R**
    - Makes parcel-wise gene expression averages for each ROI in the 200/400 7/17 Network parcellations of Scheaffer et al (2018).  
    - Also provides donor-specific maps to check that any expression effects are stable at the individual level.  
    - Since donor tissue samples vary spatially across cortex, we check that no individual donor is dominating signal in any large number of parcels.  
    - BONUS: also computes averaged expression for desikan ROIs.  
    - Will create dscalar plots for:.  
        (1) averaged PVALB/SST expression.  
        (2) number of donors with data in each parcel.  
        (3) percentage of samples in a parcel from the max donor (measure of donor over-representation).  
        (4) average sst/pvalb parcelwise plots for each donor.  
    - At the end, we also perform a control analysis to ensure that the the SST/PVALB relationship is not driven by subject sampling bias or batch effects.  
        > the 6 AHBA donors were split into 2 groups of 3. Average LH parcelwise expression of each gene was calculated for each group.  
        > then SST expression from group 1 was correlated to PVALB expression in group 2, and vice versa. Process was repeated for each possible group split.  


**06a_preprocess_lake_data.R**
    - Preprocesses raw single cell UMI visual and frontal cortex data from Lake et al. 2018.  
    - Applies basic preprocessing steps using the Seurat package.  
        (1) identification of cell and gene expression outliers.  
        (2) log-normalization.  
        (3) scaling and regression of number of detected genes per cell.  
    - Create gene name dictionary to match AHBA genes to sn-DropSeq genes.  
    - Write *seurat_processed.Rdata data objects for later reading.  


**06b_cibersortx_prep.R**
    - For visual and frontal cortex sn-DropSeq data, subset to those genes that are present in both.  
    - Transform data out of log-space (format expected by CIBERSORTx).  
    - To reduce collinearity among gene signature matrix, collapse cell subtypes into overarching categorie (e.g. In6a and In6b cell become just PVALB).  
    - Write single-cell expression matrix, and AHBA mixture files.  
    - (not shown) Polygenic deconvolution of cell type fractions using (https://cibersortx.stanford.edu/).  
        >> at the time of publication, binaries for CIBERSORTx software had not been released for integration into analysis scripts.  


**06c_cibersortx.R**
    - Process cell deconvolution output from CIBERSORTx.  
    - Plot signature matrices, and correlation plot of single-cell signature matrices for visual and frontal cortex cells.  
    - Create parcel-wise *dscalar.nii maps of cell-fractions.  
    - Correlate spatial maps of cell fractions to AHBA SST and PVALB.  


**07_get_rsfa_into_csv.R**
    - Read RSFA data from each subject and compile into a csv file.  


**08_UKB_rsfc.R**
    - Compile RSFC data, calc average connectivity for each edge.  
    - Create *dlabel.nii file showing cortical overlap between limbic striatal/default thalamus; and som/mot striatal/thalamic.  
    - Calculate SST and PVALB expression in limbic/som/mot cortical parcels;.  


**09_combine_ukb_data.R**
    - Combines UKB resting state data with phenotype info.  
    - remove subjects with any missing data in a parcel.  
    - filters subjects with run-wise motion greater than 0.20.  
    - use DOB and scan date to calc age at scan.  
    - fill in missing BMI vals.  
    - combine with lesion volume estimates.  
    - subset to white/brit subs, remove heterozygosity/invSNR outliers.  
    - write subject list:.  
        >> then execute "12a_genetics_pipe.bash" to remove cryptically related individuals and calcualate genetic PCs.  
        >> subset MRI/PHENO dataframe, combine genetic PCS, and write for later processing.  


**10_cluster_rsfa_and_write_GCTA_files.R**
    - scale quantitative covariates, residualize RSFA values for each parcel.  
    - hierarchical clustering into 7 clusters.  
    - create averaged parcel-wise RSFA average map, RSFA cluster dlabel.nii map.  
    - write csv files of phenotype, quantitative covariate, and categorical covariate data for GCTA.  
    - read AHBA and cell-fraction data for each parcel, then correlate to parcel-wise RSFA.  
    - the rest of the script produces a bunch of plots.  


**11_genetics_pipe.bash**
    - subset UKB genetic data to usable imaging subjects.  
    - basic SNP preprocessing; GENO=0.02; mind=0.1; hwe=1e-6; MAF=0.5.  
    - remove crypic relatedness and calculate eigevectors.  
    - parallel preprocessing was done with MAF=0.01. results did not change with different MAF thresholds.  


**12_make_snp.lists.R**  
    - identify SNPs occuring near coding regions of each gene  
    - read *bim SNP data form genetic preprocessing  
    - query biomaRt to retrieve information about gene locations  
    - read eQTL SNPs from NIH GTEx and CommonMind consortia  
    - cross-reference UKB genotypes SNPs with the gene dictionary, write RSID lists (i.e. for top 500 PVALB genes)  


**13_snp_herit.py**  
    - use SNP lists to make partitioned GRMs  
    - calculate overall and partitioned heritability estimates for each cluster and parcel  


**14_process_hsq_files.R**  
    - read/plot cluster-wise and parcel-wise heritability estimates (overall and partitioned)  
    - Wald test of significance for partitioned heritability estimates  
    - correlate partitioned heritabilities with parcel-wise AHBA expression, and cell-fraction maps  
    - *dscalar.nii plots of PVALB/SST partitioned heritability  


**15_prep_for_magma.bash**  
    - Run gene-based MAGMA for the Ripke SCZ GWAS  


**16_prep_for_scz_magma.R**  
    - Partitioned gene enrichment analysis of SCZ-GWAS MAGMA results  
    - Create MAGMA plots in Figure 6  


**17_scz_prs.R**  
    - remove the MHC region from the UKB genotype data (except for top p-val SNP)  
    - run PRSice to calculate SCZ polygenic risk for UKB subjects  
    - run regression of SCZ PRS to cluster-wise and parcel-wise RSFA  
    - correlate parcel-wise AHBA and cell-fraction estimates to the parcelwise map of SCZ-PRS predictions  


**18_ahba_rnaseq.R**  
    - stand-alone SST/PVALB cortical analysis of regional TPM values in the two bi-hemispheric AHBA subjects  


**19_bernard_primate_ctx.R**  
    - stand-alone SST/PVALB cortical analysis of gene expression in NIH Blueprint (Bernard et al) macaques  


**20_brainspan_dev_ctx.R**  
    - Correlation of SST and PVALB in cortex across developmental age bins using data from Brainspan Atlas.  







