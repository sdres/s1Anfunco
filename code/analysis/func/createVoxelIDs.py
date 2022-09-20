import nibabel as nb
import numpy as np

dataFile = '/home/sebastian/Desktop/testRelativeResolution/patchFlatten/sub-06_ses-001_task-stim_run-001_T1w.nii'

dataArr = nb.load(dataFile).get_fdata()

dataArr.shape
dataArrFlat = dataArr.flatten()
dataArrFlat.size
newArr = np.arange(dataArrFlat.size)
dataArrFlat.size
newArr = newArr.reshape(dataArr.shape)

newArr


header = nb.load(dataFile).header
affine = nb.load(dataFile).affine
newNii = nb.Nifti1Image(newArr, header=header,affine=affine)
nb.save(newNii,'/home/sebastian/Desktop/testRelativeResolution/patchFlatten/sub-06_ses-001_task-stim_run-001_T1w_voxelIDs.nii')
