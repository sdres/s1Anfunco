'''

Doing motion correction and regestering multiple runs across sessions.

'''

import ants
import os
import glob
from nipype.interfaces import afni
import nibabel as nb
import numpy as np
import subprocess
from IPython.display import clear_output
import nipype.interfaces.fsl as fsl
import itertools
import pandas as pd

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'
ANTSPATH = '/Users/sebastiandresbach/ANTs/install/bin'

def computeT1w(nulledFile, notnulledFile):

    '''

    Takes nulled and notnulled files as input and computes T1w image
    in EPI space. Returns array instead of saving a file to allow different
    naming conventions.

    '''

    # Load nulled motion corrected timeseries
    nulledNii = nb.load(nulledFile)
    nulledData = nulledNii.get_fdata()

    # Load notnulled motion corrected timeseries
    notnulledNii = nb.load(notnulledFile)
    notnulledData = notnulledNii.get_fdata()

    # Concatenate nulled and notnulled timeseries
    combined = np.concatenate((notnulledData,nulledData), axis=3)

    # Compute std deviation
    stdDev = np.std(combined, axis=3)
    #Compute mean
    mean = np.mean(combined, axis=3)
    # Compute variation
    cvar = stdDev/mean
    # Take inverse
    cvarInv = 1/cvar

    return cvarInv

subs = ['sub-06']

for sub in subs:
    ses = 'ses-01'

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    if not os.path.exists(funcDir):
        os.makedirs(funcDir)
        print("Func directory is created")

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    # make folder to dump motion traces
    motionDir = f'{funcDir}/motionParameters'
    if not os.path.exists(motionDir):
        os.makedirs(motionDir)
        print("Motion directory is created")

    refBase = os.path.basename(runs[0]).rsplit('.', 2)[0][:-4]


    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        # make folder to dump motion traces for the run
        os.system(f'mkdir {funcDir}/motionParameters/{base}')


        for start, modality in enumerate(['notnulled', 'nulled']):
            print(modality)
            # Load timeseries containing nulled and notnulled
            nii = nb.load(run)
            # get header and affine
            header = nii.header
            affine = nii.affine
            # Load data as array
            dataComplete = nii.get_fdata()

            # separate nulled and notnulled data
            data = dataComplete[:,:,:,start:-2:2] # Start is defined by "enumerate" above. 0 for notnulled, 1 for nulled. Here, I also get rid of the noise maps
            # make new nii and save
            img = nb.Nifti1Image(data, header=header, affine=affine)
            nb.save(img, f'{funcDir}/{base}_{modality}.nii')


            if 'rest' in run:
                # make reference image
                reference = np.mean(data[:,:,:,4:6],axis=-1)
                # and save it
                img = nb.Nifti1Image(reference, header=header, affine=affine)
                nb.save(img, f'{funcDir}/{base}_{modality}_reference.nii')

            # define mask and reference images in 'antspy-style'
            fixed = ants.image_read(f'{funcDir}/{refBase}_{modality}_reference.nii')
            mask = ants.get_mask(fixed, cleanup=2)
            ants.image_write(mask, f'{funcDir}/{base}_{modality}_moma.nii')


            # Load data in antsPy style
            ts = ants.image_read(f'{funcDir}/{base}_{modality}.nii')

            # Perform motion correction
            corrected = ants.motion_correction(ts, fixed = fixed, mask = mask)
            ants.image_write(corrected['motion_corrected'], f'{funcDir}/{base}_{modality}_moco.nii')

            # save transformation matrix for later
            for vol, matrix in enumerate(corrected['motion_parameters']):
                mat = matrix[0]
                os.system(f"cp {mat} {funcDir}/motionParameters/{base}/{base}_{modality}_vol{vol:03d}.mat")

        # =========================================================================
        # Compute T1w image in EPI space within run

        t1w = computeT1w(f'{funcDir}/{base}_nulled_moco.nii', f'{funcDir}/{base}_notnulled_moco.nii')
        header = nb.load(f'{funcDir}/{base}_nulled_moco.nii').header
        affine = nb.load(f'{funcDir}/{base}_nulled_moco.nii').affine
        img = nb.Nifti1Image(t1w, header=header, affine=affine)
        nb.save(img, f'{funcDir}/{base}_T1w.nii')


    # =========================================================================
    # Compute T1w image in EPI space across runs

    for i, run in enumerate(runs):
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        for modality in ['notnulled', 'nulled']:
            data = nb.load(f'{funcDir}/{base}_{modality}_moco.nii').get_fdata()

            if i == 0 and modality == 'notnulled':

                combined = data.copy()

                header = nb.load(f'{funcDir}/{base}_{modality}_moco.nii').header
                affine = nb.load(f'{funcDir}/{base}_{modality}_moco.nii').affine

            else:

                combined = np.concatenate((combined, data), axis = 3)

    stdDev = np.std(combined, axis=3)
    mean = np.mean(combined, axis=3)
    cvar = stdDev/mean
    cvarInv = 1 / cvar

    img = nb.Nifti1Image(cvarInv, header=header, affine=affine)
    nb.save(img, f'{funcDir}/{sub}_{ses}_T1w.nii')

    t1w = ants.image_read(f'{funcDir}/{sub}_{ses}_T1w.nii')
    t1w_n4 = ants.n4_bias_field_correction(t1w)
    ants.image_write(t1w_n4, f'{funcDir}/{sub}_{ses}_T1w_N4corrected.nii')
