import os
import numpy as np
from scipy.ndimage import morphology
from nibabel import load, save, Nifti1Image



# Load data
nii = load('/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/seg_rim_polished_perimeter_chunk.nii.gz')
basename = nii.get_filename().split(os.extsep, 1)[0]
dirname = os.path.dirname(nii.get_filename())
data = np.asarray(nii.dataobj)


data = data >= 1

data = morphology.binary_dilation(data, iterations=20)

out = Nifti1Image(data, header=nii.header, affine=nii.affine)
save(out, basename + "_dilate.nii.gz")
