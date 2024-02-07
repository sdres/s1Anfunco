"""Script to extract laminar resting state activity from digit ROIs"""

import nibabel as nb
import numpy as np
import subprocess
import os
import glob


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-07']

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest'

    # Load digit ROIs
    roiFile = f'{funcDir}/rois/{sub}_BOLD_allRois.nii.gz'
    roiNii = nb.load(roiFile)  # Load nifti
    roiData = roiNii.get_fdata()  # Load data as array

    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    layers = np.unique(depthData)[1:]

    for modality in ['VASO', 'BOLD']:
        restVolumes = sorted(glob.glob(f'{regRestDir}/vol_*_reg.nii.gz'))

        for volume in restVolumes:
            # Load rest volume
            restVolNii = nb.load(volume)  # Load nifti
            restVolData = restVolNii.get_fdata()  # Load data as array

            val = np.mean(np.multiply(restVolData, csfBrain))

            csfTimeCourse.append(val)

