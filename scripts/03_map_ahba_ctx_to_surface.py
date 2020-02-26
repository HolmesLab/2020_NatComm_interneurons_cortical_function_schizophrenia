#!/bin/python

import os
import scipy.io as sio
import nibabel as nib
import numpy as np
import pandas as pd
from nibabel import cifti2 as ci
from subprocess import Popen, PIPE

# this script will use native space AHBA sample coordinates (x,y,z) to map onto
# a 32k midthickness file that is spatially aligned with native cortical
# geometry (vertices are aligned across individuals)

# transform native voxel coordinates to RAS scanner space
def aff_conv(i, j, k, aff_matrix, aff_vector):
    return aff_matrix.dot([i, j, k]) + aff_vector

# Set up directories/paths
data_path = '/gpfs/milgram/project/holmes/kma52/sst_pvalb_nn/data'
ahba_path = os.path.join(data_path, 'ahba')
fs_path   = os.path.join(data_path, 'ahba_fs')
wb_command_path = '/gpfs/milgram/apps/hpc.rhel7/software/ConnectomeWorkbench/1.3.2/bin_rh_linux64/wb_command'


# Read sample information - AHBA
ahba_samples = pd.read_csv(os.path.join(data_path, 'ahba/ahba_sampleInfo.csv'))
RAS_df       = ahba_samples
RAS_df['R']  = np.nan
RAS_df['A']  = np.nan
RAS_df['S']  = np.nan
RAS_df['vertex'] = 0
RAS_df['mm_to_surf'] = np.nan


# perform sample to surface mapping separately for each AHBA donor
donor_arr = ['9861', '10021', '12876', '14380', '15496', '15697']
for donor_num in donor_arr:
    print(donor_num)
    
    # read AHBA donor Freesurfer nifti file
    donor_t1 = os.path.join(data_path, 'ahba_fs/donor' + donor_num + '/mri/orig/001.nii')
    img      = nib.load(donor_t1)


    # identify cortical samples (get rid of subcortex/cerebellum/etc)
    donor_samples = ahba_samples.loc[(ahba_samples['brain'] == int(donor_num))]
    donor_idxs    = ahba_samples.index[(ahba_samples['brain'] == int(donor_num))]
    donor_ctx_samples = ahba_samples.loc[
        (ahba_samples['brain'] == int(donor_num)) & (ahba_samples['top_level'] == 'CTX')]
    donor_ctx_idxs = ahba_samples.index[
        (ahba_samples['brain'] == int(donor_num)) & (ahba_samples['top_level'] == 'CTX')]


    # affine matrix information needed to convert i,j,k (AHBA voxel locations) to RAS needed for wb_command vertex mapping
    aff_matrix = img.affine[:3, :3]
    aff_vector = img.affine[:3, 3]
    T1_center  = (np.array(img.shape) - 1) / 2.

    # Create midthickness file (avg of white and pial)
    donor_surf_dir = os.path.join(fs_path, 'donor' + str(donor_num), 'surf')

    # convert voxel index to RAS coordinates (using affine matrix)
    print('Converting IJK to RAS')
    rowct = 0
    for row in donor_idxs:
        print(str(rowct) + '/' + str(len(donor_idxs)))
        cur_sample = ahba_samples.iloc[row]
        i = int(cur_sample['mri_voxel_x'])
        j = int(cur_sample['mri_voxel_y'])
        k = int(cur_sample['mri_voxel_z'])
        RAS = aff_conv(i, j, k, aff_matrix, aff_vector)
        RAS_df.loc[row, 'R'] = RAS[0]
        RAS_df.loc[row, 'A'] = RAS[1]
        RAS_df.loc[row, 'S'] = RAS[2]
        rowct += 1

    print('Mapping RAS voxels to 32k midthick vertices')
    for hemi in ['left', 'right']:
        # find intersection of cortex/donor/hemi samples
        ctx_idxs   = ahba_samples.index[ahba_samples['top_level'] == 'CTX']
        donor_idxs = ahba_samples.index[ahba_samples['brain'] == int(donor_num)]
        hemi_idxs  = ahba_samples.index[ahba_samples['structure_name'].str.contains(hemi)]
        
        # get intersection of donor/ctx/hemi indices
        donor_ctx_idx       = ctx_idxs.intersection(donor_idxs)
        donor_ctx_hemi_idxs = donor_ctx_idx.intersection(hemi_idxs)

        # LH only AHBA donors will give nulls for right hemisphere loops, so skip
        if len(donor_ctx_hemi_idxs) == 0:
            continue

        # sample data for this donor/hemi/ctx
        hemi_dat = RAS_df.loc[donor_ctx_hemi_idxs]

        # identify the correct midthick file for mapping
        if hemi == 'left':
            midthick = os.path.join(donor_surf_dir, 'lh.fsLR.32k.C_RAS.midthickness.native.surf.gii')
        elif hemi == 'right':
            midthick = os.path.join(donor_surf_dir, 'rh.fsLR.32k.C_RAS.midthickness.native.surf.gii')

        # use wb_command to identify the closest cortical vertex of each AHBA sample i,j,k location
        donor_RAS   = hemi_dat[['R', 'A', 'S']].loc[donor_ctx_hemi_idxs]
        voxel_file  = os.path.join(data_path,
                                  'ahba_fs/fsLR32k/donor' + str(donor_num) + '_' + hemi + '_RAS_fs32k_voxels.txt')
        vertex_file = os.path.join(data_path,
                                   'ahba_fs/fsLR32k/donor' + str(donor_num) + '_' + hemi + '_RAS_2_fs32k_vertices.txt')
        donor_RAS.to_csv(voxel_file, sep=' ', index=False, header=False)
        os.system(wb_command_path + ' -surface-closest-vertex ' + midthick + ' ' + voxel_file + ' ' + vertex_file)

        # read info from wb_command; i.e. closest vertex for each RAS voxel
        vertex_mapping = pd.read_csv(vertex_file, header=None, names=['vertex'])
        RAS_df.loc[donor_ctx_hemi_idxs, 'vertex'] = vertex_mapping.as_matrix()

        # create signed distance volume of every voxel to the surface
        dist_vol   = os.path.join(donor_surf_dir, hemi + '.signed_dist_vol.fsLR.32k.C_RAS.midthickness.native.nii.gz')
        signed_vol = wb_command_path + ' -create-signed-distance-volume ' + midthick + ' ' + donor_t1 + ' ' + dist_vol
        os.system(signed_vol)

        # use afni to extract distance values at each RAS coordinate
        # this will allow us to remove samples that are far from a given vertex location (i.e. less precise surface mapping)
        maskdump = '/gpfs/milgram/apps/hpc.rhel7/software/AFNI/2018.08.28/3dmaskdump'
        for row in donor_ctx_hemi_idxs:
            cur_row = RAS_df.iloc[row]
            process = Popen([maskdump, '-nbox', str(cur_row['R']), str(cur_row['A']), str(cur_row['S']), dist_vol],
                            stdout=PIPE, stderr=PIPE)
            stdout, stderr = process.communicate()
            dist_to_surf = float(str(stdout).split(' ')[-1].replace('\\n', '').replace("'", ''))
            print(dist_to_surf)
            RAS_df.loc[row, 'mm_to_surf'] = dist_to_surf

# write sample info with vertex mapping to csv
RAS_df.to_csv(os.path.join(data_path, 'ahba/sample_info_vertex_mapped.csv'), index=False)










