import numpy as np
import nibabel as nb
import os
import glob
import matplotlib.pyplot as plt
import subprocess
import os.path
from skimage.measure import label



digits=['D2','D3','D4']
subs = ['sub-02', 'sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-05']


# Threshold and binarize activation
for sub in subs:
    print(f'processing {sub}')
    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    mapFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'
    outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp'

    for digit in digits:
        print(digit)
        mapFile = f'{mapFolder}/{digit}vsAll_BOLD.nii'
        baseName = mapFile.split('/')[-1].split('.')[0]
        mapNii = nb.load(mapFile)
        mapData = mapNii.get_fdata()
        thr = np.max(mapData)/2
        mapData = mapData >= thr

        data = label(mapData, connectivity=1)
        labels, counts = np.unique(data, return_counts=True)
        largestCluster = np.argmax(counts[1:])+1

        tmp = data == largestCluster
        img = nb.Nifti1Image(tmp, affine = mapNii.affine, header = mapNii.header)
        nb.save(img, f'{outFolder}/{baseName}_largestCluster_bin.nii.gz')

# use UVD filter in LayNii to propagate ROI across cortical depth if voxel in GM
for sub in subs:

    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    mapFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'
    outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp'

    for digit in digits:
        mapFile = f'{mapFolder}/{digit}vsAll_BOLD.nii'
        baseName = mapFile.split('/')[-1].split('.')[0]

        command = '/home/sebastian/git/laynii/LN2_UVD_FILTER '
        command += f'-values {outFolder}/{baseName}_largestCluster_bin.nii.gz '
        command += f'-coord_uv {rimFolder}/seg_rim_polished_UV_coordinates.nii.gz '
        command += f'-coord_d {rimFolder}/seg_rim_polished_metric_equivol.nii.gz '
        command += f'-domain {rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz '
        command += f'-radius 0.45 '
        command += f'-height 2 '
        command += f'-max'

        subprocess.run(command,shell=True)


# sometimes there are floating bits. Therefore, we can select largest cluster of mask
for sub in subs:
    print(f'processing {sub}')
    outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp'

    for digit in digits:
        print(digit)
        file = f'{outFolder}/{digit}vsAll_BOLD_largestCluster_bin_UVD_max_filter.nii.gz'
        baseName = mapFile.split('/')[-1].split('.')[0]
        nii = nb.load(file)
        data = nii.get_fdata()


        clusterData = label(data, connectivity=1)
        labels, counts = np.unique(clusterData, return_counts=True)
        largestCluster = np.argmax(counts[1:])+1

        tmp = clusterData == largestCluster
        img = nb.Nifti1Image(tmp, affine = mapNii.affine, header = mapNii.header)
        nb.save(img, f'{outFolder}/{baseName}_largestCluster.nii.gz')


# make nii containing all rois in a binarized version. This will then be inflated and registered to the functional data. In this way we can mask the functional data o the portions that we are interested in and heavily reduce data load.
for sub in subs:
    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    mapFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'
    outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/tmp'

    allRoisNii = nb.load(f'{outFolder}/{digit}vsAll_BOLD_largestCluster_bin_UVD_max_filter.nii.gz')
    affine = allRoisNii.affine
    header = allRoisNii.header
    allRoisData = allRoisNii.get_fdata()
    allRoisData =  np.zeros(allRoisData.shape)


    for i, digit in enumerate(digits, start=1):
        roiFile = f'{outFolder}/{digit}vsAll_BOLD_largestCluster_bin_UVD_max_filter.nii.gz'
        roiNii = nb.load(roiFile)
        roiData = roiNii.get_fdata()
        tmp = roiData*i

        allRoisData = np.add(allRoisData,tmp)

        # remove voxels that overlap between ROIs
        allRoisData[allRoisData > i] = 0


    img = nb.Nifti1Image(allRoisData, affine = affine, header = header)
    nb.save(img, f'{outFolder}/{sub}_allRois.nii.gz')

    allRoisData[allRoisData > 0] = 1
    img = nb.Nifti1Image(allRoisData, affine = affine, header = header)
    nb.save(img, f'{outFolder}/{sub}_allRois_bin.nii.gz')

# for sub in subs:
#     rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
#     outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'
#     for digit in digits:
#         for modality in ['BOLD', 'VASO']:
#             mapFile = f'{outFolder}/{digit}vsAll_{modality}_flat_1000x1000_voronoi_smooth.nii.gz'
#             mapNii = nb.load(mapFile)
#             mapData = mapNii.get_fdata()
#
#             maxData = np.max(mapData)
#             thr = maxData/2
#
#             mask = np.zeros(mapData.shape)
#
#             for i in range(mapData.shape[0]):
#                 for j in range(mapData.shape[1]):
#                     if np.any(mapData[i,j,:] >= thr):
#                         mask[i,j,:] = 1
#
#
#             img = nb.Nifti1Image(mask, header = mapNii.header, affine=mapNii.affine)
#             nb.save(img, f'{outFolder}/{digit}_{modality}_roi.nii.gz')
#
#
#         command = f'fslmaths '
#         command += f'{outFolder}/{digit}_BOLD_roi.nii.gz '
#         command += f'-add {outFolder}/{digit}_VASO_roi.nii.gz '
#         command += f'-bin '
#         command += f'{outFolder}/{digit}_roi.nii.gz'
#         subprocess.run(command, shell=True)
#
#
# subList = []
# roiList = []
# stimList = []
# modalityList = []
# layerList = []
# valList = []
#
# for sub in subs:
#     rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
#     outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'
#
#     layersFile = f'{rimFolder}/seg_rim_polished_layers_equivol_flat_1000x1000_voronoi.nii.gz'
#     layersNii = nb.load(layersFile)
#     layersData = layersNii.get_fdata()
#     layers = np.unique(layersData)[1:]
#
#     for readROI in digits:
#
#         roiFile = f'{outFolder}/{readROI}_VASO_roi.nii.gz'
#         roiNii = nb.load(roiFile)
#         roiData = roiNii.get_fdata()
#
#         layersRoiData = np.multiply(roiData,layersData)
#
#         for stim in digits:
#
#
#             for modality in ['BOLD', 'VASO']:
#                 mapFile = f'{outFolder}/{stim}vsAll_{modality}_flat_1000x1000_voronoi.nii'
#                 mapNii = nb.load(mapFile)
#                 mapData = mapNii.get_fdata()
#
#                 for i in layers:  # Compute bin averages
#                     layerRoi = layersRoiData == i
#                     mask_mean = np.mean(mapData[layerRoi])
#
#                     valList.append(mask_mean)
#                     subList.append(sub)
#                     roiList.append(readROI)
#                     stimList.append(stim)
#                     modalityList.append(modality)
#                     layerList.append(i)
#
#
# import pandas as pd
# import seaborn as sns
# data = pd.DataFrame({'subject': subList, 'readRoi': roiList, 'stimDigit':stimList,'modality':modalityList, 'layer': layerList, 'data':valList})
#
#
# g = sns.FacetGrid(data, col="readRoi",  row="modality", hue='stimDigit',sharey=False)
# g.map(sns.lineplot, "layer", "data")
# g.add_legend()
#
# roiFile = f'/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/VASO_roi.nii.gz'
# roiNii = nb.load(roiFile)
# roiData = roiNii.get_fdata()
# roiIDX = roiData == 1
#
# layersFile = f'/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/seg_rim_polished_layers_equivol_flat_1000x1000_voronoi.nii.gz'
# layersNii = nb.load(layersFile)
# layersData = layersNii.get_fdata()
# layersRoiData = np.multiply(roiData,layersData)
#
# layers = np.unique(layersData)[1:]
#
# fig, ax = plt.subplots()
# for modality in ['BOLD', 'VASO']:
#     data = []
#     mapFile = f'/home/sebastian/Desktop/patchFlatten/sub-02_done/upsampled/registeredFunc/D2vsAll_{modality}_flat_1000x1000_voronoi.nii'
#     mapNii = nb.load(mapFile)
#     mapData = mapNii.get_fdata()
#     roiData = mapData[roiIDX]
#
#     for i in layers:  # Compute bin averages
#         layerRoi = layersData == i
#         mask_mean = np.mean(mapData[layerRoi])
#         data.append(mask_mean)
#     plt.plot(layers,data,label=modality)
#
# plt.show()
#
#
#
#
#
#
#
#
#
#
#
#         # test
