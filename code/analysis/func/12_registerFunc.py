"""

Register resting state data to anatomical images


"""
import nibabel as nb
import numpy as np
import subprocess
import os
import glob

FUNCDATADIR = '/Volumes/extData/S1ANFUNCO/derivatives'
ANATDATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives/patchFlatten'

# Set subjects to work on
subs = ['sub-07']

for sub in subs:

    # Make dir to temporarily dump registered volumes
    tmpDir = f'{ANATDATADIR}/{sub}_done/upsampled/registeredFunc/regTmp'
    if not os.path.exists(tmpDir):
        os.makedirs(tmpDir)
        print("Tmp directory is created")

    # Find subject-specific rest-run
    restRun = glob.glob(f'{FUNCDATADIR}/{sub}/func/{sub}_ses-00*_task-restBA3b_run-001_VASO_masked.nii.gz')[0]

    # Load rest run
    restRun = nb.load(restRun)  # Load nifti
    header = restRun.header  # Get header
    affine = restRun.affine  # Get affine
    restData = restRun.get_fdata()  # Load data as array

    # Load anatomical data
    anatFile = f'{ANATDATADIR}/{sub}_done/upsampled/anat_scaled.nii.gz'
    anatNii = nb.load(anatFile)  # Load nifti
    anatData = anatNii.get_fdata()  # Load data as array
    anatHead = anatNii.header  # Get header
    anatAff = anatNii.affine  # Get affine

    # Loop over volumes of resting-state data
    for vol in range(restData.shape[-1]):
        # Get data of current volume
        tmp = restData[:, :, :, vol]
        # Save as invididual file
        img = nb.Nifti1Image(tmp, header=header, affine=affine)
        nb.save(img, f'{tmpDir}/vol_{vol}.nii.gz')

        # Register volume to anatomical image
        command = f'antsApplyTransforms '
        command += f'--interpolation BSpline[5] '
        command += f'-d 3 -i {tmpDir}/vol_{vol}.nii.gz '
        command += f'-r {anatFile} '
        command += f'-t {ANATDATADIR}/{sub}_done/upsampled/registeredFunc/registered1_1Warp.nii.gz '
        command += f'-t {ANATDATADIR}/{sub}_done/upsampled/registeredFunc/registered1_0GenericAffine.mat '
        command += f'-o {tmpDir}/vol_{vol}_reg.nii.gz'
        subprocess.run(command, shell=True)

        # Mask with ROI
        command = f'fslmaths {tmpDir}/vol_{vol}_reg.nii.gz -mul {ANATDATADIR}/{sub}_done/upsampled/tmp/allRois_bin.nii.gz {tmpDir}/vol_{vol}_reg.nii.gz'
        subprocess.run(command, shell=True)

    # Get all volumes
    volumes = sorted(glob.glob(f'{tmpDir}/vol_*_reg.nii.gz'))
    # Get shape of anatomical data (our new reference space)
    shape = anatData.shape

    # Make new data in reference spoace but number of volumes according to func
    new = np.empty((shape[0], shape[1], shape[2], len(volumes)), dtype="float16")
    # Fill new data with registered volumes
    for i, vol in enumerate(volumes):
        tmp = nb.load(vol).get_fdata()
        new[:, :, :, i] = tmp

    # Save as new file
    img = nb.Nifti1Image(new, header=anatHead, affine=anatAff)
    nb.save(img, f'{tmpDir}/func_reg.nii.gz')


    # subprocess.run('ImageMath 4 test5.nii.gz TimeSeriesAssemble 1 0 /home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/regTmp/vol_*_reg.nii.gz')
