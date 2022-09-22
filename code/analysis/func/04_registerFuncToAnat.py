'''

Registering statistical maps to anatomical data

'''

import glob
import os
import subprocess

subs = ['sub-06']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'



for sub in subs:



    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatDir = f'{ROOT}/derivatives/{sub}/anat'
    fixed = glob.glob(f'{anatDir}/{sub}_*_highres-mp2rage_average_uni_N4corrected_trunc.nii.gz')[0]
    moving = glob.glob(f'{funcDir}/sub-06_*_T1w.nii')[-1]

    command = 'antsRegistration '
    command += f'--verbose 1 '
    command += f'--dimensionality 3 '
    command += f'--float 0 '
    command += f'--collapse-output-transforms 1 '
    command += f'--interpolation BSpline[5] '
    command += f'--output [registered1_,registered1_Warped.nii,1] '
    command += f'--use-histogram-matching 0 '
    command += f'--winsorize-image-intensities [0.005,0.995] '
    command += f'--initial-moving-transform {ROOT}/derivatives/manualSteps/{sub}/objective_matrix.txt '
    command += f'--transform SyN[0.1,3,0] '
    command += f'--metric CC[{fixed}, {moving},1,2] '
    command += f'--convergence [60x10,1e-6,10] '
    command += f'--shrink-factors 2x1 '
    command += f'--smoothing-sigmas 1x0vox '
    command += f'-x {anatDir}/sub-06_registrationMask.nii.gz'

    subprocess.run(command,shell=True)

    command = 'antsApplyTransforms '
    command += f'--interpolation BSpline[5] '
    command += f'-d 3 -i {moving} '
    command += f'-r {fixed} '
    command += f'-t registered1_1Warp.nii.gz '
    command += f'-t registered1_0GenericAffine.mat '
    command += f'-o {moving.split(".")[0]}_registered-anat.nii'

    subprocess.run(command,shell=True)
