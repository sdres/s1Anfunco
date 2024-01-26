"""
Save nuisiance regressors so that they can be used in a FSL GLM

"""

import glob
import pandas as pd
import os
import numpy as np
from scipy import interpolate
import nibabel as nb

subs = ['sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-12']

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

modalities = {'nulled': 'VASO', 'notnulled': "BOLD"}

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'
        for start, modality in enumerate(['notnulled', 'nulled']):

            # Load motion regressors
            motionRegs = pd.read_csv(f'{motionDir}/{base}_{modality}_motionRegressors.csv')

            # ====================================================================================================
            # Upsample motion regressors

            motionRegsUps = pd.DataFrame()
            for col in motionRegs.columns:
                tmp = motionRegs[col].to_numpy()
                x = np.arange(0, tmp.shape[0])
                interp = interpolate.interp1d(x, tmp, kind='linear', bounds_error=None, fill_value='extrapolate')
                xUp = np.arange(0, tmp.shape[0], 0.5)
                new = interp(xUp)
                motionRegsUps[col] = new

            # Load motion outliers
            motionOuts = np.loadtxt(f'{motionDir}/{base}_{modality}_motionOutliers.txt')

            # ====================================================================================================
            # Upsample motion outliers

            for arrNr in range(motionOuts.shape[-1]):
                tmp = motionOuts[:, arrNr]
                x = np.arange(0, tmp.shape[0])
                interp = interpolate.interp1d(x, tmp, kind='nearest', bounds_error=None, fill_value= 'extrapolate')
                xUp = np.arange(0, tmp.shape[0], 0.5)
                new = interp(xUp)

                motionRegsUps[f'scrub {arrNr}'] = new

            # Duplicate first timepoint for VASO
            if modality == 'nulled':
                # Duplicate first row
                tmp = motionRegsUps.iloc[0:1, :]
                motionRegsUps = pd.concat((tmp,motionRegsUps))

            # ================================================================================================
            # Cut dataframe to actual number of volumes

            nii = nb.load(f'{funcDir}/{base}_{modalities[modality]}.nii.gz')  # Load data
            header = nii.header  # Get header
            nVols = header['dim'][4]  # Get nr of volumes from header
            motionRegsUps = motionRegsUps.iloc[:nVols]  # Limit dataframe

            # ================================================================================================
            # Save dataframe
            motionRegsUps.to_csv(f'{motionDir}/{base}_{modalities[modality]}_nuisance.txt',
                                 index=False,
                                 sep='\t',
                                 header=None
                                 )
