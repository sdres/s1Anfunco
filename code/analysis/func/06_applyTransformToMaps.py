"""

Registering statistical maps to anatomical data

"""

import glob
import subprocess

subs = ['sub-05', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-07']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'  # Location of anatomical data
    regFolder = f'{anatDir}/registrationFiles'

    for digit in ['D2', 'D3', 'D4']:
        for modality in ['BOLD', 'VASO']:
            for contrast in ['VsRest', 'VsAll']:
                statMap = f'{funcDir}/statMaps/{sub}_{digit}{contrast}_{modality}.nii.gz'
                print(f'Registering {sub}_{digit}{contrast}_{modality}')

                fixed = f'{anatDir}/anat_scaled.nii.gz'
                moving = statMap

                command = 'antsApplyTransforms '
                command += f'--interpolation BSpline[5] '
                command += f'-d 3 '
                command += f'-i {moving} '
                command += f'-r {fixed} '
                command += f'-t {regFolder}/registered1_1Warp.nii.gz '
                command += f'-t {regFolder}/registered1_0GenericAffine.mat '
                command += f'-o {moving.split(".")[0]}_registered.nii'

                subprocess.run(command, shell=True)

# ======================================================
# Register tSNR

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'  # Location of anatomical data
    regFolder = f'{anatDir}/registrationFiles'

    for modality in ['BOLD', 'VASO']:

        fixed = f'{anatDir}/anat_scaled.nii.gz'
        moving = glob.glob(f'{funcDir}/{sub}_ses-0*_task-restBA3b_run-001_{modality}_tSNR.nii.gz')[0]
        print(f'Registering {moving}')

        command = 'antsApplyTransforms '
        command += f'--interpolation BSpline[5] '
        command += f'-d 3 '
        command += f'-i {moving} '
        command += f'-r {fixed} '
        command += f'-t {regFolder}/registered1_1Warp.nii.gz '
        command += f'-t {regFolder}/registered1_0GenericAffine.mat '
        command += f'-o {moving.split(".")[0]}_registered.nii'

        subprocess.run(command, shell=True)


# ======================================================
# Register brain mask

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'  # Location of anatomical data
    regFolder = f'{anatDir}/registrationFiles'

    fixed = f'{anatDir}/anat_scaled.nii.gz'
    moving = glob.glob(f'{funcDir}/{sub}_ses-0*_task-restBA3b_run-001_nulled_moma.nii.gz')[0]
    print(f'Registering {moving}')

    command = 'antsApplyTransforms '
    command += f'--interpolation NearestNeighbor '
    command += f'-d 3 '
    command += f'-i {moving} '
    command += f'-r {fixed} '
    command += f'-t {regFolder}/registered1_1Warp.nii.gz '
    command += f'-t {regFolder}/registered1_0GenericAffine.mat '
    command += f'-o {moving.split(".")[0]}_registered.nii'

    subprocess.run(command, shell=True)
