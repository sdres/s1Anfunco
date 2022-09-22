import subprocess
import glob
import os
import nibabel as nb
import numpy as np
import re


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

AFNIPATH = '/Users/sebastiandresbach/abin'
TR = 1.9295

for sub in ['sub-06']:
    # for ses in ['ses-001', 'ses-002']:
    for ses in ['ses-01']:

        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-*run-00*_cbv.nii.gz'))

        outFolder = f'{ROOT}/derivatives/{sub}/func'

        for run in runs:
            base = os.path.basename(run).rsplit('.', 2)[0][:-4]
            print(f'Processing run {base}')

            for start, modality in enumerate(['notnulled', 'nulled']):

                command = f'{AFNIPATH}/3dUpsample '
                command += f'-overwrite '
                command += f'-datum short '
                command += f'-prefix {outFolder}/{base}_{modality}_intemp.nii '
                command += f'-n 2 '
                command += f'-input {outFolder}/{base}_{modality}_moco.nii'

                subprocess.call(command, shell=True)

                # fix TR in header
                subprocess.call(
                    f'3drefit -TR {TR} '
                    + f'{outFolder}'
                    + f'/{base}_{modality}_intemp.nii',
                    shell=True
                    )

                # =====================================================================
                # Duplicate first nulled timepoint to match timing between cbv and bold
                # =====================================================================

                if modality == 'nulled':
                    nii = nb.load(
                        f'{outFolder}'
                        + f'/{base}_{modality}_intemp.nii'
                        )

                    # Load data
                    data = nii.get_fdata()  # Get data
                    header = nii.header  # Get header
                    affine = nii.affine  # Get affine

                    # Make new array
                    newData = np.zeros(data.shape)

                    for i in range(data.shape[-1]-1):
                        if i == 0:
                            newData[:,:,:,i] = data[:,:,:,i]
                        else:
                            newData[:,:,:,i] = data[:,:,:,i-1]

                    # Save data
                    img = nb.Nifti1Image(newData, header=header, affine=affine)
                    nb.save(img, f'{outFolder}'
                        + f'/{base}_{modality}_intemp.nii'
                        )


                # ==========================================================================
                # BOLD-correction
                # ==========================================================================

                nulledFile = f'{outFolder}/{base}_nulled_intemp.nii'
                notnulledFile = f'{outFolder}/{base}_notnulled_intemp.nii'

                # Load data
                nii1 = nb.load(nulledFile).get_fdata()  # Load cbv data
                nii2 = nb.load(notnulledFile).get_fdata()  # Load BOLD data

                # nii1 = nii1 + 1
                # nii2 += 1
                #
                # nii1[nii1 == 0]

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
                new = np.divide(nii1[:,:,:,:maxTP], nii2[:,:,:,:maxTP])

                # Clip range to -1.5 and 1.5. Values should be between 0 and 1 anyway.
                new[new > 1.5] = 1.5
                new[new < -1.5] = -1.5

                # Save bold-corrected VASO image
                img = nb.Nifti1Image(new, header=header, affine=affine)

                nb.save(
                    img, f'{outFolder}'
                    + f'/{base}_VASO_LN.nii'
                    )



            subprocess.run(f'fslmaths {outFolder}/{base}_VASO_LN.nii -mul 100 {outFolder}/{base}_VASO.nii.gz -odt short', shell=True)
            subprocess.run(f'fslmaths {outFolder}/{base}_notnulled_intemp.nii  -mul 1 {outFolder}/{base}_BOLD.nii.gz -odt short', shell=True)

            # calculate quality measures
            for modality in ['BOLD', 'VASO']:
                subprocess.run(f'LN_SKEW -input {outFolder}/{base}_{modality}.nii.gz', shell=True)
