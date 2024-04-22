"""

Upsampling registered statistical maps to match segmentation and layering

"""

import os
import glob
import subprocess

# Set data path
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-12']

for sub in subs:
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim*run-00*_cbv.nii.gz'))
    funcDir = f'{ROOT}/derivatives/{sub}/func/statMaps'

    # Define output folder
    outFolder = f'{funcDir}/upsampled'

    # Create folder if not exists
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        print("Output directory is created")

    for run in runs:

        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        digits = ['D2', 'D3', 'D4']
        if sub == 'sub-12':
            digits.append('D1')
            digits.append('D5')

        for digit in digits:
            for modality in ['BOLD', 'VASO']:
                for contrast in ['VsRest', 'VsAll']:
                    map = f'{funcDir}/{base}_{modality}_{digit}{contrast}_registered.nii'
                    outBase = os.path.basename(map).split('.')[0]

                    # Upsamling
                    command = f'c3d '
                    command += f'{map} '
                    command += f'-resample 300x300x300% '
                    command += f'-interpolation Cubic '
                    command += f'-o {outFolder}/{outBase}_upsampled.nii'

                    subprocess.run(command, shell=True)

