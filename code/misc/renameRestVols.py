"""Rename rest volumes"""
import os

import nibabel as nb
import numpy as np
import glob
import subprocess

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-14']
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-16', 'sub-17', 'sub-18']


for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest'

    for modality in ['VASO', 'BOLD']:
        restVolumes = sorted(glob.glob(f'{regRestDir}/{sub}_ses-0*_{modality}_vol*_registered.nii.gz'))

        for vol in restVolumes:
            parts = vol.split("_")
            volNr = parts[-2]

            if len(volNr) >= 6:
                continue

            if len(volNr) == 4:
                print(volNr)
                volNr = f'vol00{volNr[-1]}'
                print(volNr)

            if len(volNr) == 5:
                print(volNr)
                volNr = f'vol0{volNr[-2:]}'
                print(volNr)

            newName = ''
            for part in parts[:-2]:
                newName += f'{part}_'
            newName += f'{volNr}_registered.nii.gz'

            print(vol)
            print(newName)
            command = f'mv {vol} {newName}'

            subprocess.run(command, shell=True)
