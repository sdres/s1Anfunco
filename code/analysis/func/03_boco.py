"""

Performs temporal upsampling and BOLD-correction.
Also computes some QA measures.

"""

import subprocess
import glob
import os
import nibabel as nb
import numpy as np

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

TR = 1.9295

subs = ['sub-12']

for sub in subs:

    # Look for all runs in all sessions
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-*/func/{sub}_ses-*_task-*run-00*_cbv.nii.gz'))

    # Look for sessions of this participant
    sessions = []
    for run in runs:
        for i in range(1, 4):
            if f'ses-0{i}' in run:
                ses = f'ses-0{i}'
                sessions.append(ses)

    sessions = set(sessions)  # Remove duplicates

    for ses in sessions:
        # Look for runs of current session
        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-*run-00*_cbv.nii.gz'))
        # Set output folder
        outFolder = f'{ROOT}/derivatives/{sub}/func'

        # Loop over runs
        for run in runs:
            # Set base name of run
            base = os.path.basename(run).rsplit('.', 2)[0][:-4]
            print(f'Processing run {base}')

            for start, modality in enumerate(['notnulled', 'nulled']):

                # =============================================================================
                # Temporal upsampling

                # Prepare command
                command = f'3dUpsample '
                command += f'-overwrite '
                command += f'-datum short '
                command += f'-prefix {outFolder}/{base}_{modality}_intemp.nii.gz '
                command += f'-n 2 '
                command += f'-input {outFolder}/{base}_{modality}_moco_trunc.nii.gz'
                # Run command
                subprocess.call(command, shell=True)

                # =============================================================================
                # Fix TR in header
                subprocess.call(
                    f'3drefit -TR {TR} '
                    + f'{outFolder}'
                    + f'/{base}_{modality}_intemp.nii.gz',
                    shell=True
                    )

                # =====================================================================
                # Duplicate first nulled timepoint to match timing between cbv and bold

                if modality == 'nulled':
                    nii = nb.load(f'{outFolder}/{base}_{modality}_intemp.nii.gz')

                    # Load data
                    data = nii.get_fdata()  # Get data
                    header = nii.header  # Get header
                    affine = nii.affine  # Get affine

                    # Make new array
                    newData = np.zeros(data.shape)

                    # Loop over timepoints
                    for i in range(data.shape[-1]):
                        # Just take first timepoint as it is
                        if i == 0:
                            newData[..., i] = data[..., i]
                        # Shift all other timepoints one back
                        else:
                            newData[..., i] = data[..., i-1]

                    # Save data
                    img = nb.Nifti1Image(newData, header=header, affine=affine)
                    nb.save(img, f'{outFolder}/{base}_{modality}_intemp.nii.gz')

            # ==========================================================================
            # BOLD-correction

            nulledFile = f'{outFolder}/{base}_nulled_intemp.nii.gz'
            notnulledFile = f'{outFolder}/{base}_notnulled_intemp.nii.gz'

            # Load data
            nii1 = nb.load(nulledFile).get_fdata()  # Load cbv data
            nii2 = nb.load(notnulledFile).get_fdata()  # Load BOLD data

            # Find timeseries with fewer volumes
            if nii1.shape[-1] < nii2.shape[-1]:
                maxTP = nii1.shape[-1]
            elif nii1.shape[-1] > nii2.shape[-1]:
                maxTP = nii2.shape[-1]
            else:
                maxTP = nii1.shape[-1]-1

            header = nb.load(nulledFile).header  # Get header
            affine = nb.load(nulledFile).affine  # Get affine

            # Divide VASO by BOLD for actual BOCO
            new = np.divide(nii1[..., :maxTP], nii2[..., :maxTP])

            # Clip range to -1.5 and 1.5. Values should be between 0 and 1 anyway.
            new[new > 1.5] = 1.5
            new[new < -1.5] = -1.5

            # Save BOLD-corrected VASO image
            img = nb.Nifti1Image(new, header=header, affine=affine)
            nb.save(img, f'{outFolder}/{base}_VASO_LN.nii.gz')

            # Save data as short floats and multiply VASO by 100
            subprocess.run(f'fslmaths {outFolder}/{base}_VASO_LN.nii.gz -mul 100 {outFolder}/{base}_VASO.nii.gz -odt short', shell=True)
            subprocess.run(f'fslmaths {outFolder}/{base}_notnulled_intemp.nii.gz  -mul 1 {outFolder}/{base}_BOLD.nii.gz -odt short', shell=True)

            # Calculate quality measures
            for modality in ['BOLD', 'VASO']:
                subprocess.run(f'LN_SKEW -input {outFolder}/{base}_{modality}.nii.gz', shell=True)
