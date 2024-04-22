"""Script to extract laminar resting state activity from digit ROIs"""

import nibabel as nb
import numpy as np
import glob
import os

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
# subs = ['sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-06', 'sub-09', 'sub-10', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-15', 'sub-16', 'sub-17']
# subs = ['sub-07']

subs = ['sub-05', 'sub-12']
subs = ['sub-12']

digits = ['D2', 'D3', 'D4']

# ====================================================
# Extracting from INDEXED arrays

for sub in subs:
    print(f'Processing {sub}')

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    regRestDir = f'{funcDir}/registeredRest'

    # ====================================================
    # get shape for index
    volumes = sorted(glob.glob(f'{regRestDir}/peri/*BOLD*vol*_registered.nii.gz'))
    base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]
    vol1 = nb.load(volumes[0]).get_fdata()
    idx = vol1 > 0

    # ====================================================
    # Get depth information
    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    idxLayers = depthData[idx]
    layers = np.unique(idxLayers)

    # ====================================================
    # Load digit ROIs
    roiFile = f'{funcDir}/rois/{sub}_BOLD_allRois.nii.gz'
    roiNii = nb.load(roiFile)  # Load nifti
    roiData = roiNii.get_fdata()  # Load data as array
    idxRois = roiData[idx]

    for modality in ['VASO', 'BOLD']:
        print(f'Extracting from {modality}')
        # ====================================================
        # Load previously extracted data
        volumes = sorted(glob.glob(f'{regRestDir}/peri/*{modality}*vol*_registered.nii.gz'))
        base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]
        data = np.loadtxt(f'{regRestDir}/{base}_maskedPeri_clean.csv', delimiter=',')

        # Prepare empty data
        nrVols = data.shape[0]
        nrRegions = len(layers) * len(digits)
        timecourses = np.zeros((nrVols, nrRegions))

        for i in range(nrVols):
            if i % 10 == 0:
                print(f'Extracting from volume {i}')

            # Get rest volume
            volData = data[i, :]

            offset = 0

            for j, digit in enumerate(digits, start=1):
                roiIdx = (idxRois == j).astype("bool")

                for k, layer in enumerate(layers):
                    layerIdx = (idxLayers == layer).astype('bool')

                    tmp = roiIdx * layerIdx

                    val = np.mean(volData[tmp])

                    timecourses[i, (offset+k)] = val

                offset = offset + len(layers)

        np.savetxt(f'results/{sub}_{modality}_roi-BOLD_layerTimecourses_clean.csv', timecourses, delimiter=',')


# ====================================================
# Extracting from entire volumes

#
# for sub in subs:
#     print(f'Processing {sub}')
#     funcDir = f'{ROOT}/derivatives/{sub}/func'
#     anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
#
#     # Make dir to temporarily dump registered volumes
#     regRestDir = f'{funcDir}/registeredRest'
#
#     depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
#     depthNii = nb.load(depthFile)
#     depthData = depthNii.get_fdata()
#     layers = np.unique(depthData)[1:]
#
#     # for roiModality in ['VASO', 'BOLD']:
#     for roiModality in ['BOLD']:
#         # Load digit ROIs
#         roiFile = f'{funcDir}/rois/{sub}_{roiModality}_allRois.nii.gz'
#         roiNii = nb.load(roiFile)  # Load nifti
#         roiData = roiNii.get_fdata()  # Load data as array
#
#         for modality in ['VASO', 'BOLD']:
#             print(f'Extracting from {modality}')
#             restVolumes = sorted(glob.glob(f'{regRestDir}/{sub}_ses-0*_{modality}_vol*_registered.nii.gz'))
#             nrVols = len(restVolumes)
#             print(f'Found {nrVols} volumes')
#
#             nrRegions = len(layers) * len(digits)
#
#             timecourses = np.zeros((nrVols, nrRegions))
#             # np.savetxt(f'results/{sub}_{modality}_layerTimecoursesTest.csv', timecourses, delimiter=',')
#
#             for i, volume in enumerate(restVolumes):
#                 if i % 10 == 0:
#                     print(f'Extracting from volume {i}')
#
#                 # Load rest volume
#                 restVolNii = nb.load(volume)  # Load nifti
#                 restVolData = restVolNii.get_fdata()  # Load data as array
#                 offset = 0
#
#                 for j, digit in enumerate(digits, start=1):
#                     roiIdx = (roiData == j).astype("bool")
#
#                     for k, layer in enumerate(layers):
#                         layerIdx = (depthData == layer).astype('bool')
#
#                         tmp = roiIdx * layerIdx
#                         val = np.mean(restVolData[tmp])
#
#                         timecourses[i, (offset+k)] = val
#
#                     offset = offset + len(layers)
#
#             np.savetxt(f'results/{sub}_{modality}_roi-{roiModality}_layerTimecourses.csv', timecourses, delimiter=',')
