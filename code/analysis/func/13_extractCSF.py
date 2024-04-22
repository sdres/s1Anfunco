"""Script to extract CSF as nuosance regressor for resting state analysis"""

import nibabel as nb
import numpy as np
import subprocess
import os
import glob


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-14']

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest'

    # # Load registered brain mask data
    maskFile = f'{anatFolder}/{sub}_csfPlusPeri.nii.gz'
    maskNii = nb.load(maskFile)  # Load nifti
    maskData = maskNii.get_fdata()  # Load data as array

    periFile = f'{anatFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    periNii = nb.load(periFile)
    periData = periNii.get_fdata()
    perihead = periNii.header
    periAff = periNii.affine

    periData = np.where(periData == 1, 1, 0)

    csf = np.subtract(maskData, periData)
    new = nb.Nifti1Image(csf.astype('int'), affine=periAff, header=perihead)
    nb.save(new, f'{anatFolder}/csfTest.nii.gz')

    for modality in ['VASO', 'BOLD']:
        print(f'Processing {modality}')
        restVolumes = sorted(glob.glob(f'{regRestDir}/periCsf/*{modality}*vol*_registered.nii.gz'))
        nrVols = len(restVolumes)

        timecourse = np.zeros(nrVols)

        for i, volume in enumerate(restVolumes):
            if i % 10 == 0:
                print(f'Extracting from volume {i}')

            # Load rest volume
            restVolNii = nb.load(volume)  # Load nifti
            restVolData = restVolNii.get_fdata()  # Load data as array

            tmp = (restVolData * csf).astype('bool')

            val = np.mean(restVolData[tmp])

            timecourse[i] = val
            np.savetxt(f'results/{sub}_{modality}_roi-CSF_timecourse.csv', timecourse, delimiter=',')


# data = np.loadtxt(f'results/{sub}_{modality}_roi-CSF_timecourse.csv', delimiter=',')


# data = np.loadtxt(f'results/sub-14_BOLD_roi-BOLD_layerTimecourses.csv', delimiter=',')
