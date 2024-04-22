"""Make discrete distance bins for figure"""

import nibabel as nb
import numpy as np
import subprocess

root = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'
digits = ['D2', 'D3', 'D4']

for digit in digits:
    distFile = f'{root}/sub-07/func/rois/sub-07_{digit}_maxpoint_BOLD_geodistance.nii.gz'
    nii = nb.load(distFile)
    distData = nii.get_fdata()
    maxDist = int(np.max(distData))

    # Set distance bin based on maximum distance
    distances = np.arange(0, 14 + 1, 2)

    for i, dist in enumerate(distances[:-1], start=1):
        distData[(distData >= (distances[i-1]+0.1)) & (distData <= distances[i])] = i

    distData[distData > distances[-1]] = 0

    new = nb.Nifti1Image(distData, header=nii.header, affine=nii.affine)
    nb.save(new, f'{root}/sub-07/func/rois/sub-07_{digit}_maxpoint_BOLD_geodistance_bin.nii.gz')

    command = f'fslmaths {root}/sub-07/func/rois/sub-07_{digit}_maxpoint_BOLD_geodistance_bin.nii.gz '
    command += f'-add {root}/sub-07/anat/upsampled/gm_figure.nii.gz '
    command += f'{root}/sub-07/func/rois/sub-07_{digit}_segGeo.nii.gz'
    subprocess.run(command, shell=True)

# ================================================================
# Get rid of ROIs

for digit in digits:
    # Load distance file with segmentation

    distFile = f'{root}/sub-07/func/rois/sub-07_{digit}_segGeo.nii.gz'
    distNii = nb.load(distFile)
    distData = distNii.get_fdata()

    roiFile = f'{root}/sub-07/func/rois/sub-07_BOLD_allRois_bin.nii.gz'
    roiNii = nb.load(roiFile)
    roiData = roiNii.get_fdata()

    newData = np.where(roiData == 1, 1, distData)
    new = nb.Nifti1Image(newData, header=distNii.header, affine=distNii.affine)
    nb.save(new, f'{root}/sub-07/func/rois/sub-07_{digit}_segGeo_noROIs.nii.gz')


# ========================================================
# Check for subjects that have remaining 0-2 mm bin voxels

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-07']

for sub in subs:

    roiFile = f'{root}/{sub}/func/rois/{sub}_BOLD_allRois_bin.nii.gz'
    roiNii = nb.load(roiFile)
    roiData = roiNii.get_fdata()

    for digit in digits:
        distFile = f'{root}/{sub}/func/rois/{sub}_{digit}_maxpoint_BOLD_geodistance.nii.gz'
        nii = nb.load(distFile)
        distData = nii.get_fdata()
        maxDist = int(np.max(distData))

        # Set distance bin based on maximum distance
        distances = np.arange(0, 14 + 1, 2)

        for i, dist in enumerate(distances[:-1], start=1):
            distData[(distData >= (distances[i-1]+0.1)) & (distData <= distances[i])] = i

        distData[distData > distances[-1]] = 0
        newData = np.where(roiData == 1, 0, distData)

        firstBin = np.where(newData == 1, 1, 0)
        nrVox = np.sum(firstBin)
        if nrVox > 0:
            print(f'{sub} {digit}')

            if sub == 'sub-07':
                tmp = nb.Nifti1Image(firstBin, header=distNii.header, affine=distNii.affine)
                nb.save(tmp, f'{root}/sub-07/func/rois/sub-07_{digit}_segGeo_noROIs_firstBin.nii.gz')

# ========================================================
# Plot remaining 0-2 mm bin voxels



hexes = {'#003f5c': ['#003f5c', '#8699aa', '#ffffff'],
         '#374c80': ['#374c80', '#9ba1be', '#ffffff'],
         '#7a5195': ['#7a5195', '#bca5c9', '#ffffff'],
         '#bc5090': ['#bc5090', '#e1a9c6', '#ffffff'],
         '#ef5675': ['#ef5675', '#ffafb7', '#ffffff'],
         '#ff764a': ['#ff764a', '#ffbca2', '#ffffff'],
         '#ffa600': ['#ffa600', '#ffd291', '#ffffff']
         }

for hex in hexes:
    h = hex.lstrip('#')
    print('RGB =', tuple(int(h[i:i+2], 16) for i in (0, 2, 4)))


# Check if max point is in ROI

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-07']
digits = ['D2', 'D3', 'D4']

for sub in subs:
    roisFile = f'{root}/{sub}/func/rois/{sub}_BOLD_allRois.nii.gz'
    roisData = nb.load(roisFile).get_fdata()

    for i, digit in enumerate(digits, start=1):
        # Load max poin
        maxFile = f'{root}/{sub}/func/rois/{sub}_{digit}_maxpoint_BOLD.nii.gz'
        maxData = nb.load(maxFile).get_fdata()
        coords = np.unravel_index(np.argmax(maxData), maxData.shape)

        isIncluded = (roisData[coords[0], coords[1], coords[2]] == i)
        if not isIncluded:
            print(f'{sub} {digit}')

