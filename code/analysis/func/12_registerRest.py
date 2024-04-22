"""

Register functional resting state data to anatomical images


"""

import nibabel as nb
import numpy as np
import subprocess
import os
import glob

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-12']

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    regFolder = f'{anatFolder}/registrationFiles'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest/peri'
    if not os.path.exists(regRestDir):
        os.makedirs(regRestDir)
        print("Directory for registered rest is created")
    #
    # # Load anatomical data
    anatFile = f'{anatFolder}/anat_scaled.nii.gz'

    # # Get perimeter file
    periFile = f'{anatFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    periNii = nb.load(periFile)
    periData = periNii.get_fdata()

    perNoBorderData = np.where(periData == 1, 1, 0)
    img = nb.Nifti1Image(perNoBorderData, header=periNii.header, affine=periNii.affine)
    nb.save(img, f'{anatFolder}/seg_rim_polished_perimeter_chunk_bin.nii.gz')

    for modality in ['VASO', 'BOLD']:
        # Find subject-specific rest-run
        restRun = glob.glob(f'{funcDir}/{sub}_ses-0*_task-restBA3b_run-001_{modality}.nii.gz')[0]

        base = os.path.basename(restRun).rsplit('.', 2)[0]
        print(base)

        # Load rest run
        restRun = nb.load(restRun)  # Load nifti
        header = restRun.header  # Get header
        affine = restRun.affine  # Get affine
        restData = restRun.get_fdata()  # Load data as array

        # Loop over volumes of resting-state data
        for vol in range(restData.shape[-1]):
            # Get data of current volume
            tmp = restData[:, :, :, vol]
            # Save as individual file
            img = nb.Nifti1Image(tmp, header=header, affine=affine)
            nb.save(img, f'{regRestDir}/{base}_vol{vol:03d}.nii.gz')

            inFile = f'{regRestDir}/{base}_vol{vol:03d}.nii.gz'
            outFile = f'{regRestDir}/{base}_vol{vol:03d}_registered.nii.gz'

            if os.path.isfile(outFile):
                print(f'{outFile} exists. skipping.')
                continue

            print(f'Registering {inFile}')

            # Register volume to anatomical image
            command = f'antsApplyTransforms '
            command += f'--interpolation BSpline[5] '
            command += f'-d 3 '
            command += f'-i {inFile} '
            command += f'-r {anatFile} '
            command += f'-t {regFolder}/registered1_1Warp.nii.gz '
            command += f'-t {regFolder}/registered1_0GenericAffine.mat '
            command += f'-o {outFile}'
            subprocess.run(command, shell=True)

            # Mask with ROI
            command = f'fslmaths {outFile} '
            command += f'-mul {anatFolder}/seg_rim_polished_perimeter_chunk_bin.nii.gz '
            # command += f'-mul {anatFolder}/{sub}_csfPlusPeri.nii.gz '
            command += f'{outFile}'

            subprocess.run(command, shell=True)
