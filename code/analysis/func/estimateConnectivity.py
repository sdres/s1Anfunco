import nibabel as nb
import numpy as np
import subprocess
import os
import glob

subs = ['sub-05']
digits = ['D2', 'D3', 'D4']


for sub in subs:

    volumes = sorted(glob.glob(f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc/regTmp/vol_*_reg.nii.gz'))

    nrVols = len(volumes)

    layerFile = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/seg_rim_polished_layers_equivol.nii.gz'
    layerNii = nb.load(layerFile)
    layerData = layerNii.get_fdata()

    idxLayers = np.unique(layerData)[1:]
    nrLayers = len(idxLayers)

    nrRegions = nrLayers*len(digits)

    timcourses = np.empty((nrVols, nrRegions))

    roisFile = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/allRois.nii.gz'
    roisNii = nb.load(roisFile)
    roisData = roisNii.get_fdata()

    for i, vol in enumerate(volumes):
        print(f'processing volume nr {i+1}')
        offset = 0
        volData = nb.load(vol).get_fdata()

        for digit in digits:
            roisFile = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/{digit}vsAll_BOLD_largestCluster_bin_UVD_max_filter.nii.gz'
            roisNii = nb.load(roisFile)
            roisData = roisNii.get_fdata()

            idxDigit = roisData == 1

            volDataRoi = volData[idxDigit]
            layerDataRoi = layerData[idxDigit]

            for k in idxLayers:
                pos = (int(k+offset))-1
                timcourses[i, pos] = np.mean(volDataRoi[layerDataRoi == k])
            offset = offset+nrLayers


np.savetxt(f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp/VASO_Rest_timecourses.csv', timcourses, delimiter=',')


from nilearn.connectome import ConnectivityMeasure
correlation_measure = ConnectivityMeasure(kind='correlation')
correlation_matrix = correlation_measure.fit_transform([timcourses])[0]

# Display the correlation matrix
from nilearn import plotting
# Mask out the major diagonal
np.fill_diagonal(correlation_matrix, 0)
import matplotlib.pyplot as plt
fig = plt.figure()
fig = plotting.plot_matrix(correlation_matrix, colorbar=True,
                     vmax=0.8, vmin=-0.8)
plt.savefig('correl.png')
