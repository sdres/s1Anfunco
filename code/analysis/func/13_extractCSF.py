"""Script to extract CSF as nuosance regressor for resting state analysis"""

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

    # # Load registered brain mask data
    maskFile = glob.glob(f'{funcDir}/{sub}_ses-0*_task-restBA3b_run-001_nulled_moma_registered.nii')[0]
    maskNii = nb.load(maskFile)  # Load nifti
    maskData = maskNii.get_fdata()  # Load data as array

    # Get segmentation file
    segFile = f'{anatFolder}/seg_rim_polished.nii.gz'
    segNii = nb.load(segFile)
    segData = segNii.get_fdata()
    csfData = np.where(segData == 1, 1, 0)

    # intersect segmentation with brain mask
    csfBrain = csfData * maskData

    for modality in ['VASO', 'BOLD']:
        restVolumes = sorted(glob.glob(f'{regRestDir}/vol_*_reg.nii.gz'))

        csfTimeCourse = []

        for volume in restVolumes:
            # Load rest volume
            restVolNii = nb.load(volume)  # Load nifti
            restVolData = restVolNii.get_fdata()  # Load data as array

            val = np.mean(np.multiply(restVolData, csfBrain))

            csfTimeCourse.append(val)

