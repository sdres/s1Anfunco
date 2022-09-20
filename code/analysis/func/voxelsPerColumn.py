import nibabel as nb
import numpy as np


colFile = '/home/sebastian/Desktop/testRelativeResolution/sub-06_rim_columns1000.nii.gz'
voxelIDsFile = '/home/sebastian/Desktop/testRelativeResolution/sub-06_ses-001_task-restBA3b_run-001_T1w_voxelIDs_scaled.nii'

colArr = nb.load(colFile).get_fdata()
voxelIDsArr = nb.load(voxelIDsFile).get_fdata()
header = nb.load(colFile).header
affine = nb.load(colFile).affine
voxVolume = np.prod(header.get_zooms())



newArr = np.zeros(colArr.shape)

for column in np.unique(colArr)[1:]:
    idx = colArr == column
    # colVolume = np.sum(idx)*voxVolume
    nrEPIVoxels = np.unique(voxelIDsArr[idx]).size
    nrColVoxels = np.sum(idx)
    metric = nrColVoxels/nrEPIVoxels
    newArr[idx] = metric


newNii = nb.Nifti1Image(newArr, header=header,affine=affine)
nb.save(newNii,'/home/sebastian/Desktop/testRelativeResolution/sub-06_voxelsPerColumn.nii')
