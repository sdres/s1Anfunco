"""Regress out motion and CSF fluctuations"""

import nibabel as nb
import numpy as np
import os
import glob
from nilearn import signal
import pandas as pd

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'

# Set subjects to work on
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

subs = ['sub-18']

TR = 1.9295

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest'

    # for modality in ['BOLD']:
    for modality in ['VASO', 'BOLD']:
        print(f'Processing {modality}')

        # ================================================
        # Load csf timecourse

        # csfData = np.loadtxt(f'results/{sub}_{modality}_roi-CSF_timecourse.csv', delimiter=',')

        # ================================================
        # Load motion timecourse
        motionDataFile = glob.glob(f'{ROOT}/{sub}/func/*/{sub}*restBA3b*/*-restBA3b_run-001_{modality}_nuisance.txt')[0]

        motionData = pd.read_csv(motionDataFile, delimiter='\t', header=None)

        # ================================================
        # Merge regressors

        # motionData['csf'] = csfData

        regressors = motionData

        # ================================================
        # Load functional data and combine into timeseries

        volumes = sorted(glob.glob(f'{regRestDir}/peri/*{modality}*vol*_registered.nii.gz'))

        base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]

        vol1 = nb.load(volumes[0]).get_fdata()
        idx = vol1 > 0

        shape = vol1[idx].shape
        new = np.zeros((shape[0], len(volumes)), dtype="float16")
        #
        # for i, vol in enumerate(volumes):
        #     if i % 10 == 0:
        #         print(f'Adding volume {i}')
        #
        #     nii = nb.load(vol)
        #     tmp = np.asarray(nii.dataobj)[idx]
        #
        #     new[:, i] = tmp

        # np.savetxt(f'{regRestDir}/{base}_maskedPeri.csv', new, delimiter=',')
        new = np.loadtxt(f'{regRestDir}/{base}_maskedPeri.csv', delimiter=',')

        new = new.transpose()
        skipVols = 2
        regressorsNew = regressors.iloc[skipVols:]

        cleaned = signal.clean(new[skipVols:len(regressorsNew)+skipVols, :],
                               detrend=True,
                               standardize='psc',
                               confounds=regressorsNew,
                               high_pass=0.01,
                               t_r=TR
                               )

        np.savetxt(f'{regRestDir}/{base}_maskedPeri_clean.csv', cleaned, delimiter=',')


# ====================================================================
# Check data range

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regRestDir = f'{funcDir}/registeredRest'

    # for modality in ['BOLD']:
    for modality in ['VASO', 'BOLD']:
        print(f'Processing {modality}')
        volumes = sorted(glob.glob(f'{regRestDir}/peri/*{modality}*vol*_registered.nii.gz'))

        base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]

        cleaned = np.loadtxt(f'{regRestDir}/{base}_maskedPeri_clean.csv', delimiter=',')

        minVal = np.min(cleaned)
        maxVal = np.max(cleaned)
        meanVal = np.mean(cleaned)
        print('')
        print(f'After cleaning')
        print(f'Mean: {meanVal}')
        print(f'Range: {minVal} - {maxVal}')
        print('')


for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regStimDir = f'{funcDir}/registeredStim'

    # for modality in ['BOLD']:
    for modality in ['BOLD']:
        print(f'Processing {modality}')

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)
            cleaned = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean_immediatePsc_noconfounds.csv', delimiter=',')

            minVal = np.min(cleaned)
            maxVal = np.max(cleaned)
            meanVal = np.mean(cleaned)
            print('')
            print(f'After cleaning')
            print(f'Mean: {meanVal}')
            print(f'Range: {minVal} - {maxVal}')
            print('')