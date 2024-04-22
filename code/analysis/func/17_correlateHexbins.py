"""Control analysis for resting state"""

import nibabel as nb
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math
import subprocess
import glob
import os
from nilearn.connectome import ConnectivityMeasure

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

digits = ['D2', 'D3', 'D4']


# ==================================================================================================================
# Generate Hexbins

for sub in subs:

    # ====================================================
    # Set folders

    rimFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    roiFolder = f'{funcDir}/rois'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    regRestDir = f'{funcDir}/registeredRest'

    # ============================================
    # Get average ROI volume

    roisFile = f'{roiFolder}/{sub}_BOLD_allRois.nii.gz'
    allRoisNii = nb.load(roisFile)
    allRoisData = allRoisNii.get_fdata()

    rois = np.where(allRoisData >= 1, 1, 0)
    nrVox = np.sum(rois)
    volume = nrVox * ((0.5 / 3) ** 3)

    avgVol = volume / 3

    # ============================================
    # Get average cortical thickness

    thickFile = f'{rimFolder}/seg_rim_polished_thickness.nii.gz'
    thickNii = nb.load(thickFile)
    thickData = thickNii.get_fdata()

    perimeterFile = f'{rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    perimeterNii = nb.load(perimeterFile)
    perimeterData = perimeterNii.get_fdata()
    perimeterData = np.where(perimeterData == 1, 1, 0)

    thickData = thickData * perimeterData
    idx = thickData >= 1
    thickness = np.mean(thickData[idx])

    # ============================================
    # Get perimeter volume

    periVolume = np.sum(perimeterData) * ((0.5 / 3) ** 3)

    # ============================================
    # Calculate radius of hexbins to match their volume with that of ROIs

    # Volume of a cylinder
    # V = h *  pi * r^2

    r = math.sqrt(avgVol/(thickness * math.pi))

    # ============================================
    # Actually generate hexbins

    command = '/Users/sebastiandresbach/github/LAYNII/LN2_HEXBIN '
    command += f'-coord_uv {rimFolder}/seg_rim_polished_UV_coordinates.nii.gz '
    command += f'-radius {r}'

    subprocess.run(command, shell=True)


# ==================================================================================================================
# Extract layer timecourses from Hexbins

for sub in subs:
    # ====================================================
    # Set folders

    rimFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    roiFolder = f'{funcDir}/rois'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    regRestDir = f'{funcDir}/registeredRest'

    # ====================================================
    # get shape for index

    volumes = sorted(glob.glob(f'{regRestDir}/peri/*VASO*vol*_registered.nii.gz'))
    base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]
    vol1 = nb.load(volumes[0]).get_fdata()
    idx = vol1 > 0

    # ============================================
    # Get data from hex-bins

    hexFile = sorted(glob.glob(f'{rimFolder}/*hexbins*'))[0]
    hexNii = nb.load(hexFile)
    hexData = hexNii.get_fdata()
    idxHex = hexData[idx]

    hexLabels = np.unique(idxHex)
    # print(f'Hexlabels: {hexLabels}')

    # ====================================================
    # Get depth information

    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    idxLayers = depthData[idx]
    layers = np.unique(idxLayers)

    # Remove outermost layers
    # layers = layers[1:-1]
    print(f'Layernrs: {layers}')

    for modality in ['VASO', 'BOLD']:

        print(f'Extracting from {modality}')

        # ====================================================
        # Load previously extracted data

        volumes = sorted(glob.glob(f'{regRestDir}/peri/*{modality}*vol*_registered.nii.gz'))
        base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]

        data = np.loadtxt(f'{regRestDir}/{base}_maskedPeri_clean.csv', delimiter=',')

        # Prepare empty data
        nrVols = data.shape[0]
        nrRegions = len(layers) * len(hexLabels)
        timecourses = np.zeros((nrVols, nrRegions))

        for i in range(nrVols):
            if i % 10 == 0:
                print(f'Extracting from volume {i}')

            # Get rest volume
            volData = data[i, :]

            offset = 0

            for j, label in enumerate(hexLabels, start=1):
                roiIdx = (idxHex == int(label)).astype("bool")
                # print(f'Extracting from hex label {label}')
                for k, layer in enumerate(layers):
                    layerIdx = (idxLayers == layer).astype('bool')

                    tmp = roiIdx * layerIdx

                    nrVoxels = np.sum(tmp)
                    if i == 0:
                        print(f'extracting from {nrVoxels} for hex {label} layer {layer}')
                    val = np.mean(volData[tmp])

                    timecourses[i, (offset+k)] = val

                offset = offset + len(layers)

        np.savetxt(f'results/{sub}_{modality}_roi-BOLD_hexTimecourses_clean.csv', timecourses, delimiter=',')


# ==================================================================================================================
# Correct for NaNs

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for i, sub in enumerate(subs):
    for modality in ['VASO', 'BOLD']:
        timecourses = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_hexTimecourses_clean.csv', delimiter=',')
        nrNaNs = np.sum(np.isnan(timecourses))
        print(f'{sub}: {nrNaNs} NaNs')
        nans = np.argwhere(np.isnan(timecourses))
        nans_t = nans.transpose()
        timecourses_t = timecourses.transpose()
        print(np.unique(nans_t[1]))
        # Check if NaNs are in layer 1 or 11 for all volumes
        for val in np.unique(nans_t[1]):
            if val % 11 == 0:
                print(f'missing values in layer 1, filling with layer 2')
                timecourses_t[val, :] = timecourses_t[val + 1, :]

            elif (val + 1) % 11 == 0:
                print(f'missing values in layer 11, filling with layer 10')
                timecourses_t[val, :] = timecourses_t[val - 1, :]

            else:
                print(f'missing values in other layer')

        timecourses_new = timecourses_t.transpose()
        nrNaNs = np.sum(np.isnan(timecourses_new))
        print(f'After correction: {nrNaNs} NaNs')

        np.savetxt(f'results/{sub}_{modality}_roi-BOLD_hexTimecourses_clean.csv', timecourses_new, delimiter=',')


# ==================================================================================================================
# Create null matrices for individual participants

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for i, sub in enumerate(subs):
    for modality in ['VASO', 'BOLD']:
        # Load hexbin timecourses across layers
        timecourses = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_hexTimecourses_clean.csv', delimiter=',')

        # Correlate hexbin layers
        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([timecourses])[0]

        # =========================================================
        # Make average sub-matrix

        # Initialize empty matrix
        avgMatrix = np.zeros((11, 11))

        # Define number of layers
        layers = 11

        # Get number of sub-matrices
        nrSubMatrices = int(correlation_matrix.shape[0] / layers)

        # Set starting coordinates to loop over large matrix
        yOffSet = 0
        xOffSet = 0

        # Start counting added matrixes for averaging process
        matrixCount = 0

        # Loop over rows
        for x in range(nrSubMatrices):
            print('')
            print(f'new column (nr {x})')
            print(f'from {xOffSet} to {xOffSet + layers}')
            print('')

            # Loop over columns within current row
            for y in range(nrSubMatrices):

                print(f'from {yOffSet} to {yOffSet+layers}')
                # Select submatrix
                subMatrix = correlation_matrix[yOffSet:yOffSet+layers, xOffSet:xOffSet+layers]

                # Ignore diagonal
                if not yOffSet == xOffSet:
                    # Add submatrix to average matrix
                    avgMatrix += subMatrix
                    # Add counter for added matrix
                    matrixCount += 1

                # Print if matrix is skipped to verify that diagonal is not included
                if yOffSet == xOffSet:
                    print('skipping')

                # Make sure to end loop if last row index is reached
                if yOffSet+layers >= correlation_matrix.shape[0]:
                    break
                # Make sure to end loop if last column index is reached
                if xOffSet >= correlation_matrix.shape[0]:
                    break

                # Move y-offset to next sub-matrix
                yOffSet = yOffSet + layers

            # Jumping to new column after last row
            xOffSet = xOffSet + layers
            # Reset row index
            yOffSet = 0

        avgMatrix /= matrixCount

        # Remove outer layers
        avgMatrix = avgMatrix[1:-1, 1:-1]

        # Save null matrix
        np.savetxt(f'results/{sub}_{modality}_nullMatrix.csv', avgMatrix, delimiter=',')

# ==================================================================================================================
# Normalise null-matrix of individual participants


subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for i, sub in enumerate(subs):
    for modality in ['VASO', 'BOLD']:
        # Load subject specific average hexmatrix
        matrix = np.loadtxt(f'results/{sub}_{modality}_nullMatrix.csv', delimiter=',')

        # Normalise matrix

        # Find upper and lower thresh
        minNullNorm, maxNullNorm = np.percentile(matrix, (1, 99))

        # Clip to extremes
        matrix[matrix < minNullNorm] = minNullNorm
        matrix[matrix > maxNullNorm] = maxNullNorm

        # Normalise to 0-1 range
        matrixNorm = (matrix - minNullNorm) / (maxNullNorm - minNullNorm)

        # Add subject matrix to overall matrix
        avgMatrix += matrixNorm

        np.savetxt(f'results/{sub}_{modality}_nullMatrix_norm.csv', avgMatrix, delimiter=',')

# ==================================================================================================================
# Average null-matrix across participants

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for modality in ['VASO', 'BOLD']:
    avgMatrix = np.zeros((9, 9))
    for i, sub in enumerate(subs):
        # Load subject specific average hexmatrix
        matrix = np.loadtxt(f'results/{sub}_{modality}_nullMatrix_norm.csv', delimiter=',')

        # Add subject matrix to overall matrix
        avgMatrix += matrixNorm

    # Divide matrix values by number of participants
    avgMatrix /= len(subs)

    np.savetxt(f'results/group_{modality}_nullMatrix.csv', avgMatrix, delimiter=',')


# ==================================================================================================================
# Normalise average matrix

for modality in ['VASO', 'BOLD']:
    matrix = np.loadtxt(f'results/group_{modality}_nullMatrix.csv', delimiter=',')

    # Find upper and lower thresh
    minNullNorm, maxNullNorm = np.percentile(matrix, (1, 99))

    # Clip to extremes
    matrix[matrix < minNullNorm] = minNullNorm
    matrix[matrix > maxNullNorm] = maxNullNorm

    # Normalise to 0-1 range
    matrixNorm = (matrix - minNullNorm) / (maxNullNorm - minNullNorm)

    # Save normalised null matrix
    np.savetxt(f'results/group_{modality}_nullMatrix_normalised.csv', matrixNorm, delimiter=',')


# ============================================================================================
# Plotting

colors = {'BOLD': 'Reds_r',
          'VASO': 'Blues_r'}

plt.style.use('dark_background')

locs = [0.5, 4.5, 8.5]
labels = ['Deep', 'Middle', 'Superficial']

# ==============================
# Plot average normalised matrix
for modality in ['VASO', 'BOLD']:
    matrix = np.loadtxt(f'results/group_{modality}_nullMatrix_normalised.csv', delimiter=',')

    fig, ax = plt.subplots()
    sns.heatmap(matrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
    ax.set_yticks(locs, labels, fontsize=12)
    ax.set_xticks(locs, labels, fontsize=12)
    ax.set_title(f'Group {modality} null-matrix normalised', fontsize=16)

    plt.tight_layout()
    plt.savefig(f'results/group_nullMatrix_{modality}_norm.png', bbox_inches="tight")
    plt.close()

# ==============================
# Plot subject specific matrices

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for i, sub in enumerate(subs):
    for modality in ['VASO', 'BOLD']:

        matrix = np.loadtxt(f'results/{sub}_{modality}_nullMatrix.csv', delimiter=',')

        fig, ax = plt.subplots()
        sns.heatmap(matrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
        ax.set_yticks(locs, labels, fontsize=12)
        ax.set_xticks(locs, labels, fontsize=12)
        ax.set_title(f'{sub} {modality} null-matrix', fontsize=16)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_nullMatrix_{modality}.png', bbox_inches="tight")
        plt.close()




# subs = ['sub-12']
# timecourses.shape
# timecourses.shape[1]/11
#
# locs = np.arange(5, timecourses.shape[1], 11)
# labels = []
# for hexNr in range(1, int(timecourses.shape[1]/11)+1):
#     labels.append(f'Hex {int(hexNr):02d}\nlayers')
#
# locs = [locs[0], locs[-1]]
# labels = [labels[0], labels[-1]]
# len(labels)

locs = np.arange(5, timecourses.shape[1], 11)
labels = []
for hexNr in range(1, int(timecourses.shape[1] / 11) + 1):
    labels.append(f'Hex {int(hexNr):02d}\nlayers')

locs = [locs[0], locs[-1]]
labels = [labels[0], labels[-1]]





# ============================================================================================
# Retired Code
#
#
# ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'
# for i, sub in enumerate(subs):
#     funcDir = f'{ROOT}/{sub}/func'
#     anatFolder = f'{ROOT}/{sub}/anat/upsampled'
#     # Make dir to temporarily dump registered volumes
#     regRestDir = f'{funcDir}/registeredRest'
#     for modality in ['VASO', 'BOLD']:
#         volumes = sorted(glob.glob(f'{regRestDir}/peri/*{modality}*vol*_registered.nii.gz'))
#         base = os.path.basename(volumes[0]).rsplit('.', 2)[0][:-18]
#         timecourses = np.loadtxt(f'{regRestDir}/{base}_maskedPeri_clean.csv', delimiter=',')
#         nrNaNs = np.sum(np.isnan(timecourses))
#         print(f'{sub}: {nrNaNs} NaNs')


    # ============================================
    # Check volume of hexbins

    # hexFile = sorted(glob.glob(f'{rimFolder}/*hexbins*'))[0]
    # hexNii = nb.load(hexFile)
    # hexData = hexNii.get_fdata()
    #
    # hexLabels = np.unique(hexData)[1:]
    # hexvolumes = []
    # for label in hexLabels:
    #     tmp = np.where(hexData == label, 1, 0)
    #     nrVox = np.sum(tmp)
    #     volume = nrVox * ((0.5 / 3) ** 3)
    #     hexvolumes.append(volume)
    # meanHexVol = np.mean(hexvolumes)
    # volumesAll.append(meanHexVol)


# =========================================================
# Plot subject average null matrices

# locs = [0.5, 5.5, 10.5]
# labels = ['Deep', 'Middle', 'Superficial']
#
# # Native
# fig, ax = plt.subplots()
# sns.heatmap(avgMatrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{sub} {modality} "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/{sub}_nullMatrix_{modality}.png', bbox_inches="tight")
# plt.close()
#
# # Normalised
# fig, ax = plt.subplots()
# sns.heatmap(avgMatrixNorm, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{sub} {modality} normalised "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/{sub}_nullMatrixNorm_{modality}_after.png', bbox_inches="tight")
# plt.close()
#
# # Add subject matrix to average matrix across participants
# nullMatrixAll += avgMatrix
# nullMatrixAllNorm += avgMatrixNorm
#
# # Divide average matrix by number of participants going into it
# nullMatrixAll /= len(subs)
# nullMatrixAllNorm /= len(subs)
#
# # if modality == 'BOLD':
# #     vals = [0.04, 0.14]
# # if modality == 'VASO':
# #     vals = [0.03, 0.08]
#
# fig, ax = plt.subplots()
# # sns.heatmap(nullMatrixAll, ax=ax, annot=False, vmin=vals[0], vmax=vals[1], cbar=True, square=True, cmap=colors[modality])
# sns.heatmap(nullMatrixAll, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{modality} "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/group_nullMatrix_{modality}.png', bbox_inches="tight")
# plt.close()
#
# fig, ax = plt.subplots()
# # sns.heatmap(nullMatrixAll, ax=ax, annot=False, vmin=vals[0], vmax=vals[1], cbar=True, square=True, cmap=colors[modality])
# sns.heatmap(nullMatrixAllNorm, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{modality} "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/group_nullMatrix_{modality}_normBEFORE.png', bbox_inches="tight")
# plt.close()
#
# # =========================================================
# # Normalise average matrix AFTERWARDS
#
# # Find upper and lower thresh
# minNullNorm, maxNullNorm = np.percentile(nullMatrixAll, (1, 99))
#
# # Clip to extremes
# nullMatrixAll[nullMatrixAll < minNullNorm] = minNullNorm
# nullMatrixAll[nullMatrixAll > maxNullNorm] = maxNullNorm
#
# # Normalise to 0-1 range
# avgMatrixNorm = (nullMatrixAll - minNullNorm) / (maxNullNorm - minNullNorm)
#
# fig, ax = plt.subplots()
# # sns.heatmap(nullMatrixAll, ax=ax, annot=False, vmin=vals[0], vmax=vals[1], cbar=True, square=True, cmap=colors[modality])
# sns.heatmap(avgMatrixNorm, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{modality} "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/group_nullMatrix_{modality}_normAFTER.png', bbox_inches="tight")
# plt.close()
#
# # =========================================================
# # Normalise NORMALISED average matrix AFTERWARDS
#
# # Find upper and lower thresh
# minNullNorm, maxNullNorm = np.percentile(nullMatrixAllNorm, (1, 99))
#
# # Clip to extremes
# nullMatrixAllNorm[nullMatrixAllNorm < minNullNorm] = minNullNorm
# nullMatrixAllNorm[nullMatrixAllNorm > maxNullNorm] = maxNullNorm
#
# # Normalise to 0-1 range
# avgMatrixNorm = (nullMatrixAllNorm - minNullNorm) / (maxNullNorm - minNullNorm)
#
# fig, ax = plt.subplots()
# # sns.heatmap(nullMatrixAll, ax=ax, annot=False, vmin=vals[0], vmax=vals[1], cbar=True, square=True, cmap=colors[modality])
# sns.heatmap(avgMatrixNorm, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
# ax.set_yticks(locs, labels, fontsize=12)
# ax.set_xticks(locs, labels, fontsize=12)
# ax.set_title(f'{modality} "null-matrix"', fontsize=16)
# plt.tight_layout()
# plt.savefig(f'results/group_nullMatrix_{modality}_normBEFOREANDAFTER.png', bbox_inches="tight")
# plt.close()


# ==================================================================================================================
# Create null matrices for individual participants while normalising individual matrices

# subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
#
# for i, sub in enumerate(subs):
#     for modality in ['VASO', 'BOLD']:
#         # Load hexbin timecourses across layers
#         timecourses = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_hexTimecourses_clean.csv', delimiter=',')
#
#         # Correlate hexbin layers
#         correlation_measure = ConnectivityMeasure(kind='correlation')
#         correlation_matrix = correlation_measure.fit_transform([timecourses])[0]
#
#         # =========================================================
#         # Make average sub-matrix
#
#         # Initialize empty matrix
#         avgMatrix = np.zeros((11, 11))
#
#         # Define number of layers
#         layers = 11
#
#         # Get number of sub-matrices
#         nrSubMatrices = int(correlation_matrix.shape[0] / layers)
#
#         # Set starting coordinates to loop over large matrix
#         yOffSet = 0
#         xOffSet = 0
#
#         # Start counting added matrixes for averaging process
#         matrixCount = 0
#
#         # Loop over rows
#         for x in range(nrSubMatrices):
#             print('')
#             print(f'new column (nr {x})')
#             print(f'from {xOffSet} to {xOffSet + layers}')
#             print('')
#
#             # Loop over columns within current row
#             for y in range(nrSubMatrices):
#
#                 print(f'from {yOffSet} to {yOffSet+layers}')
#                 # Select submatrix
#                 subMatrix = correlation_matrix[yOffSet:yOffSet+layers, xOffSet:xOffSet+layers]
#
#                 # Normalise matrix
#
#                 # Find upper and lower thresh
#                 minNullNorm, maxNullNorm = np.percentile(subMatrix, (1, 99))
#
#                 # Clip to extremes
#                 subMatrix[subMatrix < minNullNorm] = minNullNorm
#                 subMatrix[subMatrix > maxNullNorm] = maxNullNorm
#
#                 # Normalise to 0-1 range
#                 subMatrix = (subMatrix - minNullNorm) / (maxNullNorm - minNullNorm)
#
#                 # Ignore diagonal
#                 if not yOffSet == xOffSet:
#                     # Add submatrix to average matrix
#                     avgMatrix += subMatrix
#                     # Add counter for added matrix
#                     matrixCount += 1
#
#                 # Print if matrix is skipped to verify that diagonal is not included
#                 if yOffSet == xOffSet:
#                     print('skipping')
#
#                 # Make sure to end loop if last row index is reached
#                 if yOffSet+layers >= correlation_matrix.shape[0]:
#                     break
#                 # Make sure to end loop if last column index is reached
#                 if xOffSet >= correlation_matrix.shape[0]:
#                     break
#
#                 # Move y-offset to next sub-matrix
#                 yOffSet = yOffSet + layers
#
#             # Jumping to new column after last row
#             xOffSet = xOffSet + layers
#             # Reset row index
#             yOffSet = 0
#
#         avgMatrix /= matrixCount
#
#         # Remove outer layers
#         avgMatrix = avgMatrix[1:-1, 1:-1]
#
#         # Save null matrix
#         np.savetxt(f'results/{sub}_{modality}_nullMatrix_normaliseIndividualMats.csv', avgMatrix, delimiter=',')