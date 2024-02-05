"""

Registering T1w image in functional EPI space to high-res MP2RAGE UNI.

"""

import os
import subprocess

subs = ['sub-05']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

for sub in subs:

    # Defining folders
    funcDir = f'{ROOT}/derivatives/{sub}/func'  # Location of functional data
    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'  # Location of anatomical data

    regFolder = f'{anatDir}/registrationFiles'  # Folder where output will be saved
    # Make output folder if it does not exist already
    if not os.path.exists(regFolder):
        os.makedirs(regFolder)

    # =========================================================================
    # Registration
    # =========================================================================

    fixed = f'{anatDir}/anat_scaled.nii.gz'
    moving = f'{funcDir}/{sub}_ses-01_T1w_trunc.nii.gz'
    if sub == 'sub-14':
        moving = f'{funcDir}/{sub}_ses-03_T1w_trunc.nii.gz'
    if sub == 'sub-05':
        moving = f'{funcDir}/{sub}_ses-02_T1w_trunc.nii.gz'

    # Set up ants command
    command = 'antsRegistration '
    command += f'--verbose 1 '
    command += f'--dimensionality 3 '
    command += f'--float 0 '
    command += f'--collapse-output-transforms 1 '
    command += f'--interpolation BSpline[5] '
    command += f'--output [{regFolder}/registered1_,{regFolder}/registered1_Warped.nii,1] '
    command += f'--use-histogram-matching 0 '
    command += f'--winsorize-image-intensities [0.005,0.995] '
    command += f'--initial-moving-transform {anatDir}/registeredFunc/objective_matrix.txt '
    command += f'--transform SyN[0.1,3,0] '
    command += f'--metric CC[{fixed}, {moving},1,2] '
    command += f'--convergence [60x10,1e-6,10] '
    command += f'--shrink-factors 2x1 '
    command += f'--smoothing-sigmas 1x0vox '
    command += f'-x {anatDir}/{sub}_registrationMask.nii.gz'

    # Run command
    subprocess.run(command, shell=True)

    # Prepare command to apply transform and check quality
    command = 'antsApplyTransforms '
    command += f'--interpolation BSpline[5] '
    command += f'-d 3 -i {moving} '
    command += f'-r {fixed} '
    command += f'-t {regFolder}/registered1_1Warp.nii.gz '
    command += f'-t {regFolder}/registered1_0GenericAffine.mat '
    command += f'-o {moving.split(".")[0]}_registered.nii'
    # Run command
    subprocess.run(command, shell=True)
