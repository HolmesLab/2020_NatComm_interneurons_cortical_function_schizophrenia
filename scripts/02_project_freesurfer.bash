#!/bin/bash

# this script will project fsLR_32k maps to individual AHBA freesurfer space

module load FreeSurfer/6.0.0

# Set up directories
project_dir=/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn
AHBA_fs_dir=${project_dir}/data/ahba_fs
wb_command=/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command
wb_shortcut=/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_shortcuts

# Create midthickness, align to fs_LR_32k sphere, resample label files
for donor in 9861 10021 12876 14380 15496 15697; do
  echo ${donor}

  # FS space vol/surf data is slightly offset from native scanner coordinate space
  # This calcualtes affine transform, which is important for later vol-to-surf mapping of scanner space RAS coordinates
  MatrixX=$(mri_info ${AHBA_fs_dir}/donor${donor}/mri/brain.finalsurfs.mgz | grep "c_r" | cut -d "=" -f 5 | sed s/" "/""/g)
  MatrixY=$(mri_info ${AHBA_fs_dir}/donor${donor}/mri/brain.finalsurfs.mgz | grep "c_a" | cut -d "=" -f 5 | sed s/" "/""/g)
  MatrixZ=$(mri_info ${AHBA_fs_dir}/donor${donor}/mri/brain.finalsurfs.mgz | grep "c_s" | cut -d "=" -f 5 | sed s/" "/""/g)
  echo "1 0 0 ""$MatrixX" >  ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat
  echo "0 1 0 ""$MatrixY" >> ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat
  echo "0 0 1 ""$MatrixZ" >> ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat
  echo "0 0 0 1" >> ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat

  # -------------
  for hemi in lh rh; do
    if [ "${hemi}" == 'lh' ]; then
      structure=CORTEX_LEFT
      hemi_2=L
    elif [ "${hemi}" == 'rh' ]; then
      structure=CORTEX_RIGHT
      hemi_2=R
    fi

    # define paths to native freesurfer white/pial/sphere files
    native_white=${AHBA_fs_dir}/donor${donor}/surf/${hemi}.white
    native_pial=${AHBA_fs_dir}/donor${donor}/surf/${hemi}.pial
    native_sphere=${AHBA_fs_dir}/donor${donor}/surf/${hemi}.sphere
    native_spherereg=${AHBA_fs_dir}/donor${donor}/surf/${hemi}.sphere.reg
    native_fsLR32k_midthick=${AHBA_fs_dir}/donor${donor}/surf/${hemi}.fsLR.32k.mesh.midthickness.surf.gii

    # convert freesurfer formated surface to gifti
    mris_convert ${native_white} ${native_white}.native.surf.gii
    mris_convert ${native_pial} ${native_pial}.native.surf.gii

    # convert aparc.annot to label.gii format
    mris_convert --annot ${AHBA_fs_dir}/donor${donor}/label/${hemi}.aparc.annot ${native_white} ${AHBA_fs_dir}/donor${donor}/label/${hemi}.aparc.label.gii

    # resample the mesh to fsavg6, but keep native cortical geometry
    # -------------
    mesh_dir=/gpfs/milgram/project/holmes/HOLMES_UKB/external/Pipelines/global/templates/standard_mesh_atlases/resample_fsaverage
    new_sphere=${mesh_dir}/fs_LR-deformed_to-fsaverage.${hemi_2}.sphere.32k_fs_LR.surf.gii
    wb_shortcuts -freesurfer-resample-prep \
      ${native_white} ${native_pial} \
      ${native_spherereg} ${new_sphere} \
      ${AHBA_fs_dir}/donor${donor}/surf/${donor}.${hemi}.midthickness.surf.gii \
      ${AHBA_fs_dir}/donor${donor}/surf/${donor}.${hemi}.midthickness.32k_fs_LR.surf.gii \
      ${AHBA_fs_dir}/donor${donor}/surf/${hemi}.sphere.reg.surf.gii

    wb_command -surface-resample \
      ${AHBA_fs_dir}/donor${donor}/surf/${donor}.${hemi}.midthickness.surf.gii \
      ${AHBA_fs_dir}/donor${donor}/surf/${hemi}.sphere.reg.surf.gii \
      ${new_sphere} BARYCENTRIC \
      ${native_fsLR32k_midthick}


    # apply affine transformation
    wb_command -surface-apply-affine ${native_white}.native.surf.gii ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat ${native_white}.C_RAS.native.surf.gii
    wb_command -surface-apply-affine ${native_pial}.native.surf.gii ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat ${native_pial}.C_RAS.native.surf.gii
    wb_command -surface-apply-affine ${native_fsLR32k_midthick} ${AHBA_fs_dir}/donor${donor}/mri/c_ras.mat ${AHBA_fs_dir}/donor${donor}/surf/${hemi}.fsLR.32k.C_RAS.midthickness.native.surf.gii

    done

    # convert orig.mgz to .nii
    mri_convert ${AHBA_fs_dir}/donor${donor}/mri/orig/001.mgz ${AHBA_fs_dir}/donor${donor}/mri/orig/001.nii
done




#end
