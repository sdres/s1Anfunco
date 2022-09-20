import nibabel as nb
import numpy as np
import subprocess
import os
import glob
from scipy.ndimage import morphology


root = '/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'
subs = ['sub-05']


for sub in subs:
    # first, we have to inflate our ROIs
    # Load
    nii = nb.load(f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/{sub}_allRois_bin.nii.gz')
    basename = nii.get_filename().split(os.extsep, 1)[0]
    dirname = os.path.dirname(nii.get_filename())
    data = np.asarray(nii.dataobj)

    data = data >= 1

    data = morphology.binary_dilation(data, iterations=20)

    out = nb.Nifti1Image(data, header=nii.header, affine=nii.affine)
    nb.save(out, basename + "_dilate.nii.gz")



# register the dilated ROIs to functional data using the inverse transform
command = f'antsApplyTransforms '
command += f'--interpolation GenericLabel '
command += f'-d 3 '
command += f'-i /home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/{sub}_allRois_bin_dilate.nii.gz '
command += f'-r /media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS/derivatives/{sub}/func/{sub}_EPI_T1w_N4Corrected.nii '
command += f'-t /home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc/registered1_1InverseWarp.nii.gz '
command += f'-t [ /home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc/registered1_0GenericAffine.mat, 1] '
command += f'-o /home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/allRois_bin_dilate_toFunc.nii.gz'
subprocess.run(command, shell=True)


inFile = f'{root}/derivatives/{sub}/func/sub-05_ses-002_task-restBA3b_run-001_VASO.nii.gz'
mask = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/allRois_bin_dilate_toFunc.nii.gz'
outFile = f'{root}/derivatives/{sub}/func/sub-05_ses-002_task-restBA3b_run-001_VASO_masked.nii.gz'
# mask the functional data
command = f'fslmaths {inFile} -mul {mask} {outFile}'
subprocess.run(command, shell=True)

# now register the masked func data to our anatomical data
funcFolder = '/home/sebastian/Desktop/patchFlatten/sub-05_done/upsampled/registeredFunc'


command = 'antsApplyTransforms '
command += '--interpolation BSpline[5] '
command += '-d 4 '
command += f'-i /media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS/derivatives/{sub}/func/sub-05_ses-002_task-restBA3b_run-001_VASO_masked.nii.gz '
command += f'-r /home/sebastian/Desktop/patchFlatten/sub-05_done/upsampled/anat_scaled.nii.gz '
command += f'-t {funcFolder}/registered1_1Warp.nii.gz '
command += f'-t {funcFolder}/registered1_0GenericAffine.mat '
command += f'-o {funcFolder}/registeredRestFunc.nii.gz'

subprocess.run(command, shell=True)





for sub in subs:

    tmpDir = '/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/regTmp'
    # subprocess.run(f'mkdir {tmpDir}', shell=True)

    restRun = f'{root}/derivatives/{sub}/func/sub-02_ses-003_task-restBA3b_run-001_VASO_masked.nii.gz'
    restRun = nb.load(restRun)
    header = restRun.header
    affine = restRun.affine
    restData = restRun.get_fdata()

    anatFile = '/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/anat_scaled.nii.gz'
    anatNii = nb.load(anatFile)
    anatData = anatNii.get_fdata()
    anatHead = anatNii.header
    anatAff = anatNii.affine

    for vol in range(restData.shape[-1]):
        tmp = restData[:,:,:,vol]
        img = nb.Nifti1Image(tmp, header=header,affine=affine)
        nb.save(img, f'{tmpDir}/vol_{vol}.nii.gz')

        command = f'antsApplyTransforms '
        command += f'--interpolation BSpline[5] '
        command += f'-d 3 -i {tmpDir}/vol_{vol}.nii.gz '
        command += f'-r {anatFile} '
        command += f'-t /home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/registered1_1Warp.nii.gz '
        command += f'-t /home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/registered1_0GenericAffine.mat '
        command += f'-o {tmpDir}/vol_{vol}_reg.nii.gz'
        subprocess.run(command,shell=True)

        command = f'fslmaths {tmpDir}/vol_{vol}_reg.nii.gz -mul /home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/tmp/allRois_bin.nii.gz {tmpDir}/vol_{vol}_reg.nii.gz'
        subprocess.run(command,shell=True)

    #
    # volumes = sorted(glob.glob(f'/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/regTmp/vol_*_reg.nii.gz'))
    # shape = anatData.shape
    # new = np.empty((shape[0],shape[1],shape[2],len(volumes)), dtype="float16")
    # for i, vol in enumerate(volumes):
    #     tmp = nb.load(vol).get_fdata()
    #     new[:,:,:,i]=tmp
    # img = nb.Nifti1Image(new, header=anatHead,affine=anatAff)
    #
    # nb.save(img, '/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/regTmp/func_reg.nii.gz')
    #

    # subprocess.run('ImageMath 4 test5.nii.gz TimeSeriesAssemble 1 0 /home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/regTmp/vol_*_reg.nii.gz')
