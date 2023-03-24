"""

Doing motion correction and registering multiple runs across sessions.

"""
import ants
import os
import glob
import nibabel as nb
import numpy as np
import subprocess
from scipy import signal
import sys

sys.path.append(f'./code/misc')
from computeT1w import *

subs = ['sub-12']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    if not os.path.exists(funcDir):
        os.makedirs(funcDir)
        print("Func directory is created")

    # Look for individual runs (containing both nulled and not-nulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    for i in range(1,4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    # Make folder to dump motion traces
    motionDir = f'{funcDir}/motionParameters'
    if not os.path.exists(motionDir):
        os.makedirs(motionDir)
        print("Motion directory is created")

    # Set resting-state run as reference
    refBase = os.path.basename(runs[0]).rsplit('.', 2)[0][:-4]

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        # Make folder to dump motion traces for the run
        os.system(f'mkdir {funcDir}/motionParameters/{base}')

        for start, modality in enumerate(['notnulled', 'nulled']):

            print(modality)
            # Load timeseries containing both nulled and notnulled
            nii = nb.load(run)
            # get header and affine
            header = nii.header
            affine = nii.affine
            # Load data as array
            dataComplete = nii.get_fdata()

            # Separate nulled and not-nulled data
            # Start is defined by "enumerate" above. 0 for notnulled, 1 for nulled.
            data = dataComplete[:, :, :, start:-2:2]  # Here, I also get rid of the noise maps
            # Make new nii and save
            img = nb.Nifti1Image(data, header=header, affine=affine)
            nb.save(img, f'{funcDir}/{base}_{modality}.nii')

            if 'rest' in run:
                # Make reference image
                reference = np.mean(data[..., 4:6], axis=-1)
                # And save it
                img = nb.Nifti1Image(reference, header=header, affine=affine)
                nb.save(img, f'{funcDir}/{base}_{modality}_reference.nii')

                # Make motion mask
                print('Generating mask')
                command = '3dAutomask '
                command += f'-prefix {funcDir}/{base}_{modality}_moma.nii.gz '
                command += f'-peels 3 -dilate 2  '
                command += f'{funcDir}/{base}_{modality}_reference.nii'
                subprocess.run(command, shell=True)

            # Define mask and reference images in 'antspy-style'
            fixed = ants.image_read(f'{funcDir}/{refBase}_{modality}_reference.nii')

            # Get motion mask
            mask = ants.image_read(f'{funcDir}/{refBase}_{modality}_moma.nii.gz')

            # Load data in antsPy style
            ts = ants.image_read(f'{funcDir}/{base}_{modality}.nii')

            # Perform motion correction
            corrected = ants.motion_correction(ts, fixed=fixed, mask=mask)
            ants.image_write(corrected['motion_corrected'], f'{funcDir}/{base}_{modality}_moco.nii.gz')

            # Save transformation matrix for later
            for vol, matrix in enumerate(corrected['motion_parameters']):
                mat = matrix[0]
                os.system(f"cp {mat} {funcDir}/motionParameters/{base}/{base}_{modality}_vol{vol:03d}.mat")

        # =========================================================================
        # Compute T1w image in EPI space within run

        nulledFile = f'{funcDir}/{base}_nulled_moco.nii.gz'
        notNulledFile = f'{funcDir}/{base}_notnulled_moco.nii.gz'

        t1w = computeT1w(nulledFile, notNulledFile)
        header = nb.load(nulledFile).header
        affine = nb.load(nulledFile).affine
        img = nb.Nifti1Image(t1w, header=header, affine=affine)
        nb.save(img, f'{funcDir}/{base}_T1w.nii')

    # =========================================================================
    # Compute T1w image in EPI space across runs

    for i, run in enumerate(runs):
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        for modality in ['notnulled', 'nulled']:
            data = nb.load(f'{funcDir}/{base}_{modality}_moco.nii.gz').get_fdata()

            # For first run, initiate combined data
            if i == 0 and modality == 'notnulled':

                combined = data.copy()
                # And get header/ affine information to save data later
                header = nb.load(f'{funcDir}/{base}_{modality}_moco.nii.gz').header
                affine = nb.load(f'{funcDir}/{base}_{modality}_moco.nii.gz').affine

            # For other runs, just concatenate
            else:
                combined = np.concatenate((combined, data), axis=3)

    # De-trend before std. dev. calculation
    combinedDemean = signal.detrend(combined, axis=3, type='constant')
    combinedDetrend = signal.detrend(combinedDemean, axis=3, type='linear')
    stdDev = np.std(combinedDetrend, axis=3)

    mean = np.mean(combined, axis=3)

    cvar = np.divide(stdDev, mean)

    cvarInv = 1 / cvar

    img = nb.Nifti1Image(cvarInv, header=header, affine=affine)
    nb.save(img, f'{funcDir}/{sub}_{ses}_T1w.nii')

    t1w = ants.image_read(f'{funcDir}/{sub}_{ses}_T1w.nii')
    t1w_n4 = ants.n4_bias_field_correction(t1w)
    ants.image_write(t1w_n4, f'{funcDir}/{sub}_{ses}_T1w_N4corrected.nii')
