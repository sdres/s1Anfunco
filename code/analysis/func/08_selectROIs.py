"""

Defining ROI masks for all digits. This is done in several steps:
- Select voxels from activation map based on threshold. Here, we chose half the
maximum z-score as cut-off
- Propagate active voxels across cortical depth using LayNii's UVD filter
- Remove floating bits that might occur in the process
- Combine ROIs into one mask

"""

import matplotlib.pyplot as plt
import numpy as np
import nibabel as nb
import subprocess
import os.path
from skimage.measure import label
import pandas as pd
import seaborn as sns

# Set data path
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

SUBS = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
SUBS = ['sub-17']
DIGITS = ['D2', 'D3', 'D4']

# =============================================================================
# Threshold and binarize activation
# =============================================================================

for sub in SUBS:
    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    mapFolder = f'{funcDir}/statMaps'
    outFolder = f'{funcDir}/rois'

    # Create folder if not exists
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        print("Output directory is created")

    periFile = f'{anatDir}/seg_rim_polished_perimeter_chunk.nii.gz'
    periNii = nb.load(periFile)
    periData = periNii.get_fdata()
    periData = periData == 1

    np.unique(periData)
    for digit in DIGITS:
        for modality in ['BOLD', 'VASO']:
            for contrast in ['VsRest', 'VsAll']:

                mapFile = f'{mapFolder}/{sub}_{digit}{contrast}_{modality}_registered.nii'
                outBase = os.path.basename(mapFile).split('.')[0]

                mapNii = nb.load(mapFile)
                mapData = mapNii.get_fdata()
                thr = np.max(mapData)/3
                mapData = mapData * periData
                mapData = mapData >= thr

                data = label(mapData, connectivity=1)
                labels, counts = np.unique(data, return_counts=True)
                largestCluster = np.argmax(counts[1:])+1

                tmp = data == largestCluster
                img = nb.Nifti1Image(tmp, affine=mapNii.affine, header=mapNii.header)
                nb.save(img, f'{outFolder}/{outBase}_largestCluster_bin.nii.gz')


# =============================================================================
# Propagate ROI across cortical depth if voxel in GM
# =============================================================================

for sub in SUBS:

    anatDir = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    outFolder = f'{funcDir}/rois'
    # for modality in ['VASO', 'BOLD']:
    for modality in ['BOLD']:
        for digit in DIGITS:

            print(f'Starting with {digit}')

            mapFile = f'{outFolder}/{sub}_{digit}VsAll_{modality}_registered_largestCluster_bin.nii.gz'
            baseName = mapFile.split('/')[-1].split('.')[0]

            command = 'LN2_UVD_FILTER '
            command += f'-values {mapFile} '
            command += f'-coord_uv {anatDir}/seg_rim_polished_UV_coordinates.nii.gz '
            command += f'-coord_d {anatDir}/seg_rim_polished_metric_equivol.nii.gz '
            command += f'-domain {anatDir}/seg_rim_polished_perimeter_chunk.nii.gz '
            command += f'-radius 0.45 '
            command += f'-height 2 '
            command += f'-max'

            subprocess.run(command, shell=True)
            print(f'Done with {digit}')

# =============================================================================
# Remove floating bits
# =============================================================================

for sub in SUBS:
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    outFolder = f'{funcDir}/rois'
    # for modality in ['VASO', 'BOLD']:
    for modality in ['BOLD']:

        for digit in DIGITS:
            print(digit)

            file = f'{outFolder}/{sub}_{digit}VsAll_{modality}_registered_largestCluster_bin_UVD_max_filter.nii.gz'
            baseName = file.split('/')[-1].split('.')[0]
            nii = nb.load(file)
            data = nii.get_fdata()

            clusterData = label(data, connectivity=1)
            labels, counts = np.unique(clusterData, return_counts=True)
            largestCluster = np.argmax(counts[1:])+1

            tmp = clusterData == largestCluster
            img = nb.Nifti1Image(tmp, affine=nii.affine, header=nii.header)
            nb.save(img, f'{outFolder}/{baseName}_largestCluster.nii.gz')


# =============================================================================
# Combine ROIs within modality
# =============================================================================

for sub in SUBS:
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    outFolder = f'{funcDir}/rois'
    # for modality in ['VASO', 'BOLD']:
    for modality in ['BOLD']:
        # Create empty dataset with dimensions of target data
        file = f'{outFolder}/{sub}_D2VsAll_{modality}_registered_largestCluster_bin_UVD_max_filter_largestCluster.nii.gz'
        allRoisNii = nb.load(file)
        affine = allRoisNii.affine
        header = allRoisNii.header
        allRoisData = allRoisNii.get_fdata()
        allRoisData = np.zeros(allRoisData.shape)
        # Create array to check for overlap between ROIs
        overlap = np.zeros(allRoisData.shape)

        # Fill empty dataset with digit ROI data
        for i, digit in enumerate(DIGITS, start=1):
            roiFile = f'{outFolder}/{sub}_{digit}VsAll_{modality}_registered_largestCluster_bin_UVD_max_filter_largestCluster.nii.gz'
            roiNii = nb.load(roiFile)
            roiData = roiNii.get_fdata()
            tmp = roiData*i

            allRoisData = np.add(allRoisData, tmp)

            # Get overlap between ROIs
            tmpOverlap = np.where(allRoisData > i, 1, 0)
            overlap += tmpOverlap

            # remove voxels that overlap between ROIs
            allRoisData[allRoisData > i] = 0

        img = nb.Nifti1Image(allRoisData, affine=affine, header=header)
        nb.save(img, f'{outFolder}/{sub}_{modality}_allRois.nii.gz')

        # Save overlap
        img = nb.Nifti1Image(overlap, affine=affine, header=header)
        nb.save(img, f'{outFolder}/{sub}_{modality}_allRois_overlap.nii.gz')

        allRoisData[allRoisData > 0] = 1
        img = nb.Nifti1Image(allRoisData, affine=affine, header=header)
        nb.save(img, f'{outFolder}/{sub}_{modality}_allRois_bin.nii.gz')


# =============================================================================
# Combine ROIs across modality
# =============================================================================
#
# # make nii containing all rois in a binarized version.
# SUBS = ['sub-18']
#
# for sub in SUBS:
#     funcDir = f'{ROOT}/derivatives/{sub}/func'
#     outFolder = f'{funcDir}/rois'
#
#     # Create empty dataset with dimensions of target data
#     file = f'{outFolder}/{sub}_D2VsAll_BOLD_registered_largestCluster_bin_UVD_max_filter_largestCluster.nii.gz'
#     allRoisNii = nb.load(file)
#     affine = allRoisNii.affine
#     header = allRoisNii.header
#     allRoisData = allRoisNii.get_fdata()
#     allRoisData = np.zeros(allRoisData.shape)
#     # Create array to check for overlap between ROIs
#     roiOverlap = np.zeros(allRoisData.shape)
#
#     # Fill empty dataset with digit ROI data
#     for i, digit in enumerate(DIGITS, start=1):
#         roiFile = f'{outFolder}/{sub}_{digit}VsAll_BOLD_registered_largestCluster_bin_UVD_max_filter_largestCluster.nii.gz'
#         roiNii = nb.load(roiFile)
#         roiData = roiNii.get_fdata()
#         tmp = roiData*i
#
#         allRoisData = np.add(allRoisData, tmp)
#
#         # Get overlap between ROIs
#         tmpRoiOverlap = np.where(allRoisData > i, 1, 0)
#         roiOverlap += tmpRoiOverlap
#
#         # remove voxels that overlap between ROIs
#         allRoisData[allRoisData > i] = 0
#
#     img = nb.Nifti1Image(allRoisData, affine=affine, header=header)
#     nb.save(img, f'{outFolder}/{sub}_BOLD_allRois.nii.gz')
#
#     # Save overlap
#     img = nb.Nifti1Image(roiOverlap, affine=affine, header=header)
#     nb.save(img, f'{outFolder}/{sub}_BOLD_allRois_overlap.nii.gz')
#
#     allRoisData[allRoisData > 0] = 1
#     img = nb.Nifti1Image(allRoisData, affine=affine, header=header)
#     nb.save(img, f'{outFolder}/{sub}_BOLD_allRois_bin.nii.gz')

# =============================================================================
# Get ROI sizes (with layers)
# =============================================================================

subList = []
roiList = []
nrVoxList = []
layerList = []
roiModalityList = []

for sub in SUBS:
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    outFolder = f'{funcDir}/rois'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # for modality in ['VASO', 'BOLD']:
    for modality in ['BOLD']:

        # Load depth file
        depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
        depthNii = nb.load(depthFile)
        depthData = depthNii.get_fdata()
        layers = np.unique(depthData)[1:]

        # Load ROI file
        roisFile = f'{outFolder}/{sub}_{modality}_allRois.nii.gz'
        # roisFile = f'{outFolder}/{sub}_{digit}VsAll_{modality}_registered_largestCluster_bin_UVD_max_filter.nii.gz'
        allRoisNii = nb.load(roisFile)
        allRoisData = allRoisNii.get_fdata()

        for i, digit in enumerate(DIGITS, start=1):

            roi = np.where(allRoisData == i, 1, 0)

            for layer in layers:
                tmpLayer = np.where(depthData == layer, 1, 0)

                nrVox = np.sum(roi * tmpLayer)

                subList.append(sub)
                roiList.append(digit)
                nrVoxList.append(nrVox)
                layerList.append(layer)
                roiModalityList.append(modality)

volumeList = []
for val in nrVoxList:
    volume = val * ((0.5/3)**3)
    volumeList.append(volume)


data = pd.DataFrame({'subject': subList,
                     'layer': layerList,
                     'roi': roiList,
                     'nrVox': nrVoxList,
                     'volume': volumeList,
                     'roiModality': roiModalityList
                     })

data.to_csv(f'results/groupRoiSizesWithLayers.csv', index=False)
data = pd.read_csv(f'results/groupRoiSizesWithLayers.csv')

fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data, ax=ax, x='layer', hue='roi', y='volume', palette=paletteDigits)
# ax.set_title('ROI volume', fontsize=16)
ax.set_xticks(range(0, 11), range(1, 12), fontsize=12)
ax.set_xlabel('Layer (1 = deep)', fontsize=12)
ax.set_yticks(range(0, 15, 2), range(0, 15, 2), fontsize=12)
ax.set_ylabel('Volume [mm\u00b3]', fontsize=14)
plt.tight_layout()
plt.savefig(f'results/group_roiSizeLayers.png')
plt.show()


# =======================================================
# Without layers



subList = []
roiList = []
nrVoxList = []
roiModalityList = []

for sub in SUBS:
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    outFolder = f'{funcDir}/rois'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # for modality in ['VASO', 'BOLD']:
    for modality in ['BOLD']:

        # Load depth file
        depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
        depthNii = nb.load(depthFile)
        depthData = depthNii.get_fdata()
        layers = np.unique(depthData)[1:]

        # Load ROI file
        roisFile = f'{outFolder}/{sub}_{modality}_allRois.nii.gz'
        # roisFile = f'{outFolder}/{sub}_{digit}VsAll_{modality}_registered_largestCluster_bin_UVD_max_filter.nii.gz'
        allRoisNii = nb.load(roisFile)
        allRoisData = allRoisNii.get_fdata()

        for i, digit in enumerate(DIGITS, start=1):

            roi = np.where(allRoisData == i, 1, 0)

            nrVox = np.sum(roi)

            subList.append(sub)
            roiList.append(digit)
            nrVoxList.append(nrVox)
            roiModalityList.append(modality)

volumeList = []
for val in nrVoxList:
    volume = val * ((0.5/3)**3)
    volumeList.append(volume)

data = pd.DataFrame({'subject': subList,
                     'roi': roiList,
                     'nrVox': nrVoxList,
                     'volume': volumeList,
                     'roiModality': roiModalityList
                     })

data.to_csv(f'results/groupRoiSizes.csv', index=False)
data = pd.read_csv(f'results/groupRoiSizes.csv')


plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data, ax=ax, x='roi', y='volume', palette=paletteDigits)
ax.set_title('ROI volume', fontsize=16)
ax.set_xticks(range(3), DIGITS, fontsize=12)
ax.set_xlabel('', fontsize=12)
ax.set_yticks(range(0, 141, 20), range(0, 141, 20), fontsize=12)
ax.set_ylabel('Volume [mm\u00b3]', fontsize=14)
plt.tight_layout()
plt.savefig(f'results/group_roiSize.png')
plt.show()

# ROI sizes of participants
fig, ax = plt.subplots(figsize=(8, 5))
sns.barplot(data, ax=ax, x='subject', y='volume', hue='roi', palette=paletteDigits)
ax.set_title('ROI volume', fontsize=16)
ax.set_xticks(range(len(SUBS)), SUBS, fontsize=12)
ax.set_xlabel('', fontsize=12)
ax.set_yticks(range(0, 141, 20), range(0, 141, 20), fontsize=12)
ax.set_ylabel('Volume [mm\u00b3]', fontsize=14)
plt.tight_layout()
plt.savefig(f'results/group_roiSizePerParticipant.png')
plt.show()


for digit in DIGITS:
    tmp = data.loc[(data['roi'] == digit)]
    print(digit)
    mean = np.mean(tmp['volume'])
    print(f'Mean: {mean}')
    sd = np.std(tmp['volume'])
    print(f'Std. Dev: {sd}')
    min = 200
    max = 0
    for sub in data['subject'].unique():
        tmp = data.loc[(data['subject'] == sub) & (data['roi'] == digit)]
        mean = np.mean(tmp['volume'])
        if mean <= min:
            min = mean
        if mean >= max:
            max = mean
    print(f'lowest: {min}')
    print(f'Highest: {max}')
