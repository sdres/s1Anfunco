"""

Doing motion correction and registering multiple runs across sessions.

"""

import ants
import os
import glob
import nibabel as nb
import numpy as np
import subprocess
import sys

sys.path.append(f'./code/misc')
from computeT1w import *

subs = ['sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-07', 'sub-09', 'sub-10']
subs = ['sub-14']

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'


# ========================================================================================
# Do motion-correction within runs
print('Doing motion-correction within runs')

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/derivatives/{sub}/func'

    if not os.path.exists(funcDir):
        os.makedirs(funcDir)
        print("Func directory is created")

    # Look for individual runs (containing both nulled and not-nulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    for i in range(1, 4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    # Make folder to dump motion traces
    motionDir = f'{funcDir}/motionParameters'
    if not os.path.exists(motionDir):
        os.makedirs(motionDir)
        print("Motion directory is created")

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
            data = dataComplete[:, :, :, start::2]
            # Make new nii and save
            img = nb.Nifti1Image(data, header=header, affine=affine)
            nb.save(img, f'{funcDir}/{base}_{modality}.nii')

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
            fixed = ants.image_read(f'{funcDir}/{base}_{modality}_reference.nii')

            # Get motion mask
            mask = ants.image_read(f'{funcDir}/{base}_{modality}_moma.nii.gz')

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
            # Computing motion outliers

            command = 't '
            command += f'-i {funcDir}/{base}_{modality}.nii '
            command += f'-o {funcDir}/motionParameters/{base}/{base}_{modality}_motionOutliers.txt '
            command += f'--nomoco'
            subprocess.run(command, shell=True)

        # =========================================================================
        # Compute T1w image in EPI space within run

        nulledFile = f'{funcDir}/{base}_nulled_moco.nii.gz'
        notNulledFile = f'{funcDir}/{base}_notnulled_moco.nii.gz'

        t1w = computeT1w(nulledFile, notNulledFile)
        header = nb.load(nulledFile).header
        affine = nb.load(nulledFile).affine
        img = nb.Nifti1Image(t1w, header=header, affine=affine)
        nb.save(img, f'{funcDir}/{base}_T1w.nii')

# =============================================================================
# Register run-wise T1w images to rest
# =============================================================================

# # Set subjects to work on
# subs = ['sub-05']
# # # Set sessions to work on
# # sessions = ['ses-01', 'ses-02']
# # sessions = ['ses-05']
#
print('Registering run-wise T1w images to rest')
for sub in subs:
    print(f'Processing {sub}')
    # Look for individual STIMULATION runs (containing both nulled and not-nulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim*run-00*_cbv.nii.gz'))

    for i in range(1, 4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    outFolder = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs
    runs = sorted(glob.glob(f'{outFolder}/{sub}_{ses}_task-stim*_run-0*_T1w.nii'))

    # Set name of reference image
    refImage = (sorted(glob.glob(f'{outFolder}/{sub}_{ses}_task-restBA3b_run-001_T1w.nii')))[0]

    # Get basename of reference image
    refBase = os.path.basename(refImage).rsplit('.', 2)[0]

    # Load reference image in antsPy style
    fixed = ants.image_read(refImage)

    # Define motion mask
    mask = ants.image_read(
        sorted(glob.glob(f'{outFolder}/{sub}_{ses}_task-restBA3b_run-001_nulled_moma.nii.gz'))[0]
    )

    for run in runs:

        base = os.path.basename(run).rsplit('.', 2)[0].split('_')
        tmp = base[0]
        for subString in base[1:-1]:
            tmp = tmp + f'_{subString}'
        base = tmp

        # Define moving image
        moving = ants.image_read(run)

        # Compute transformation matrix
        mytx = ants.registration(
            fixed=fixed,
            moving=moving,
            type_of_transform='Rigid',
            mask=mask
        )

        # Apply transformation
        mywarpedimage = ants.apply_transforms(
            fixed=fixed,
            moving=moving,
            transformlist=mytx['fwdtransforms'],
            interpolator='bSpline'
        )

        # Save image
        ants.image_write(
            mywarpedimage,
            f'{outFolder}'
            + f'/{base}_T1w_registered-{refBase}.nii')

        # Get transformation name
        transform1 = (
                f'{outFolder}'
                + f'/{base}_T1w_registered-{refBase}.mat'
        )

        # Save transform for future
        os.system(f"cp {mytx['fwdtransforms'][0]} {transform1}")

# =============================================================================
# Apply between run registration
# =============================================================================
print('Apply between run registration')

# # Set subjects to work on
# subs = ['sub-05']
# # Set sessions to work on
# sessions = ['ses-01', 'ses-02']
# sessions = ['ses-04']

for sub in subs:
    print(sub)
    # Look for individual STIMULATION runs (containing both nulled and not-nulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim*run-00*_cbv.nii.gz'))

    for i in range(1, 4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    outFolder = f'{ROOT}/derivatives/{sub}/func'

    # Set name of reference image
    refImage = (sorted(glob.glob(f'{outFolder}/{sub}_{ses}_task-restBA3b_run-001_T1w.nii')))[0]

    # Get basename of reference image
    refBase = os.path.basename(refImage).rsplit('.', 2)[0].split('_')
    tmp1 = refBase[0]
    for subString in refBase[1:-1]:
        tmp1 = tmp1 + f'_{subString}'

    refBase = tmp1

    # Get header and affine of reference image to easily assign later
    refHeader = nb.load(refImage).header
    refAffine = nb.load(refImage).affine

    # Load reference image in antspy style
    fixed = ants.image_read(refImage)

    for modality in ['nulled', 'notnulled']:

        # Load motion mask of reference image in antspy style
        mask = ants.image_read(
            f'{outFolder}/{refBase}_{modality}_moma.nii.gz'
        )

        # Look for individual runs
        runs = sorted(glob.glob(f'{outFolder}/{sub}_{ses}_task-stim_run-0*_{modality}.nii'))

        # Loop over runs
        for run in runs:

            # Get the base name of the run
            runBase = os.path.basename(run).rsplit('.', 2)[0].split('_')
            tmp2 = runBase[0]
            for subString in runBase[1:-1]:
                tmp2 = tmp2 + f'_{subString}'
            runBase = tmp2

            # Load registration matrix between run and reference
            transformBetween = f'{outFolder}/{runBase}_T1w_registered-{refBase}_T1w.mat'

            # =============================================================
            # Separate run into individual volumes
            nii = nb.load(run)
            # get header and affine
            header = nii.header
            affine = nii.affine
            # Load data as array
            data = nii.get_fdata()

            # Loop over volumes
            for i in range(data.shape[-1]):
                # Overwrite volumes 0,1,2 with volumes 3,4,5
                if i <= 2:
                    vol = data[..., i + 3]
                else:
                    vol = data[..., i]

                # Save individual volumes
                img = nb.Nifti1Image(vol, header=header, affine=affine)
                nb.save(img, f'{outFolder}/{runBase}_{modality}_vol{i:03d}.nii')

            # Loop over the volumes we just created to do the correction
            for i in range(data.shape[-1]):
                # Load volume
                moving = ants.image_read(f'{outFolder}/{runBase}_{modality}_vol{i:03d}.nii')

                # Get within run transformation matrix of the volume
                transformWithin = f'{outFolder}/motionParameters' \
                                  f'/{runBase}/{runBase}_{modality}_vol{i:03d}.mat'

                # Apply transformation matrices
                mywarpedimage = ants.apply_transforms(
                    fixed=fixed,
                    moving=moving,
                    transformlist=[transformWithin, transformBetween],
                    interpolator='bSpline'
                )

                # Save warped image
                ants.image_write(mywarpedimage, f'{outFolder}/{runBase}_{modality}_vol{i:03d}_warped.nii')

            # =============================================================
            # Assemble images of run
            newData = np.zeros(data.shape)
            for i in range(data.shape[-1]):
                vol = nb.load(
                    f'{outFolder}'
                    + f'/{runBase}_{modality}_vol{i:03d}_warped.nii'
                ).get_fdata()

                newData[..., i] = vol

            img = nb.Nifti1Image(newData, header=refHeader, affine=refAffine)
            nb.save(img, f'{outFolder}/{runBase}_{modality}_moco-reg.nii')

            # Remove individual volumes
            os.system(f'rm {outFolder}/{runBase}_{modality}_vol*.nii')

    # =========================================================================
    # Compute T1w image in EPI space within run

    nulledFile = f'{outFolder}/{runBase}_nulled_moco-reg.nii'
    notNulledFile = f'{outFolder}/{runBase}_notnulled_moco-reg.nii'

    t1w = computeT1w(nulledFile, notNulledFile)
    header = nb.load(nulledFile).header
    affine = nb.load(nulledFile).affine
    img = nb.Nifti1Image(t1w, header=header, affine=affine)
    nb.save(img, f'{outFolder}/{runBase}_T1w_reg.nii')


# =========================================================================
# Compute T1w image in EPI space across runs
print('Computing T1w image in EPI space across runs')
for sub in subs:
    print(sub)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim*run-00*_cbv.nii.gz'))

    for i in range(1, 4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    outFolder = f'{ROOT}/derivatives/{sub}/func'

    for i, run in enumerate(runs):
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        for modality in ['notnulled', 'nulled']:
            data = nb.load(f'{outFolder}/{base}_{modality}_moco-reg.nii').get_fdata()

            # For first run, initiate combined data
            if i == 0 and modality == 'notnulled':

                combined = data.copy()
                # And get header/ affine information to save data later
                header = nb.load(f'{outFolder}/{base}_{modality}_moco-reg.nii').header
                affine = nb.load(f'{outFolder}/{base}_{modality}_moco-reg.nii').affine

            # For other runs, just concatenate
            else:
                combined = np.concatenate((combined, data), axis=3)

    # De-trend before std. dev. calculation
    # combinedDemean = signal.detrend(combined, axis=3, type='constant')
    # combinedDetrend = signal.detrend(combinedDemean, axis=3, type='linear')

    stdDev = np.std(combined, axis=3)

    mean = np.mean(combined, axis=3)

    cvar = np.divide(stdDev, mean)

    cvarInv = 1 / cvar

    img = nb.Nifti1Image(cvarInv, header=header, affine=affine)
    nb.save(img, f'{outFolder}/{sub}_{ses}_T1w.nii')

    # t1w = ants.image_read(f'{outFolder}/{sub}_{ses}_T1w.nii')
    # t1w_n4 = ants.n4_bias_field_correction(t1w)
    # ants.image_write(t1w_n4, f'{outFolder}/{sub}_{ses}_T1w_N4corrected.nii')
