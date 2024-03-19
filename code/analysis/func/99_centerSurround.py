"""Anna Devor inspired center surround analysis"""

import nibabel as nb
import numpy as np
import subprocess
import os
import glob
from nilearn import signal
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

subs = ['sub-16']
subs = ['sub-14', 'sub-15', 'sub-17', 'sub-18']

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12']

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'
digits = ['D2', 'D3', 'D4']

# ============================================================================
# Define center surround architecture

# for sub in subs:
#     print(f'Processing {sub}')
#     funcDir = f'{ROOT}/{sub}/func'
#     anatFolder = f'{ROOT}/{sub}/anat/upsampled'
#
#     # Load perimeter
#     perimeterFile = f'{anatFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
#     perimeterNii = nb.load(perimeterFile)
#     perimeterData = perimeterNii.get_fdata()
#     perimeterData = np.where(perimeterData == 1, 1, 0)
#
#     for digit in digits:
#         for modality in ['BOLD', 'VASO']:
#             # Load statistical map
#             statFile = f'{funcDir}/statMaps/{sub}_{digit}vsAll_{modality}_registered.nii'
#             statData = nb.load(statFile).get_fdata()
#
#             # Intersect perimeter with statistical map
#             inter = np.multiply(statData, perimeterData)
#
#             # Get index of maximum statistical voxel within gray matter
#             maxPointCoords = np.unravel_index(np.argmax(inter), statData.shape)
#
#             # Make new data and fill only maximum point
#             maxPoint = np.zeros(inter.shape)
#             maxPoint[maxPointCoords[0], maxPointCoords[1], maxPointCoords[2]] = 1
#
#             # Save maximum point
#             maxPointFile = f'{funcDir}/rois/{sub}_{digit}_maxpoint_{modality}.nii.gz'
#             ni_img = nb.Nifti1Image(maxPoint.astype('int'), affine=perimeterNii.affine, header=perimeterNii.header)
#             nb.save(ni_img, f'{maxPointFile}')
#
#             # Generate geodesic distance file
#             command = f'LN2_GEODISTANCE -domain {perimeterFile} -init {maxPointFile}'
#             subprocess.run(command, shell=True)
#

# ============================================================================
# Extract timecourses for equidistant rings


layerCompartments = {0: [1, 4],
                     1: [4, 7],
                     2: [7, 10]
                     }

subs = ['sub-14', 'sub-15']
for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'

    # get shape for index
    volumes = sorted(glob.glob(f'{funcDir}/registeredStim/peri/*BOLD*vol*_registered.nii.gz'))
    vol1 = nb.load(volumes[0]).get_fdata()
    idx = vol1 > 0

    # Get depth information
    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    idxLayers = depthData[idx]

    layers = np.unique(idxLayers)

    for modality in ['BOLD']:

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)
            if modality == 'VASO':
                if sub == 'sub-17':
                    if 'run-001' in base:
                        print('Skipping run-01 VASO for sub-17')
                        continue
            # Load perimeter timecourse data
            data = np.loadtxt(f'{funcDir}/registeredStim/{base}_maskedPeri_clean_psc.csv', delimiter=',')
            data = data.transpose()
            nrVols = data.shape[1]
            # data.shape

            for digit in digits:
                # Get maximum distance
                distFile = f'{funcDir}/rois/{sub}_{digit}_maxpoint_{modality}_geodistance.nii.gz'
                distData = nb.load(distFile).get_fdata()
                distIdx = distData[idx]
                maxDist = int(np.max(distData))

                # Set distance bin based on maximum distance
                distances = np.arange(0, maxDist+1, 2)

                # Make empty data for distance bins x timepoints
                nrRegions = len(distances)-1

                for c in layerCompartments:
                    firstLay = layerCompartments[c][0]
                    lastLay = layerCompartments[c][1]

                    mask_4 = idxLayers >= firstLay
                    mask_5 = idxLayers <= lastLay
                    mask_6 = np.logical_and(mask_4, mask_5)

                    timecourses = np.zeros((nrVols, nrRegions))
                    timecourses.shape

                    for vol in range(nrVols):
                        tmp = data[:, vol]

                        offset = 0

                        for i, distance in enumerate(distances[1:]):
                            # Mask distances
                            mask_1 = distIdx >= distances[i]
                            mask_2 = distIdx <= distance
                            mask_3 = np.logical_and(mask_1, mask_2)

                            mask = np.multiply(mask_3, mask_6)

                            val = np.mean(tmp[mask])

                            timecourses[vol, i] = val

                    np.savetxt(f'results/{base}_distFrom{digit}_layer{c+1}.csv', timecourses, delimiter=',')


# ===========================================================================
# Extract distance dependent timecourses
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
layerNames = ['Deep', 'Middle', 'Superficial']

TR = 1.9295

timePointList = []
modalityList = []
valList = []
distFromList = []
subList = []
distValList = []
runList = []
stimList = []
trialList = []
layerList = []

for sub in subs:
    print(f'Processing {sub}')

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'
    regStimDir = f'{funcDir}/registeredStim'

    for modality in ['BOLD']:
    # for modality in ['VASO', 'BOLD']:
        print(f'Extracting from {modality}')

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)
            if modality == 'VASO':
                if sub == 'sub-17':
                    if 'run-001' in base:
                        print('Skipping run-01 VASO for sub-17')
                        continue

            for i, distFromDigit in enumerate(digits):

                for c in layerCompartments:

                    # ====================================================
                    # Load previously extracted data
                    print('Load data')
                    data = np.loadtxt(f'results/{base}_distFrom{distFromDigit}_layer{c+1}.csv', delimiter=',')
                    data = data.transpose()

                    distLevels = data.shape[0]

                    # Load stim info
                    for j in range(distLevels):

                        for stimDigit in digits:

                            tmpFile = f'{ROOT}/designFiles/stimulationPatterns/{sub}/{base[:-5]}_{stimDigit}.txt'
                            tmp = pd.read_csv(tmpFile, names=['start', 'dur', 'mod'], sep=' ')

                            starts = tmp['start'].to_numpy()
                            startsTR = starts / TR
                            startsTR = np.round(startsTR).astype('int') - 15
                            ends = tmp['start'].to_numpy() + 59
                            endsTR = (ends / TR)
                            endsTR = np.round(endsTR).astype('int')

                            nrTimepoints = endsTR[0] - startsTR[0]

                            for k, (start, end) in enumerate(zip(startsTR, endsTR), start=1):
                                tmpData = data[j, start: end]
                                # tmpData = np.mean(tmpData, axis=0)

                                data.shape

                                if modality == 'VASO':
                                    tmpData = tmpData * -1

                                for h, val in enumerate(tmpData):

                                    subList.append(sub)
                                    distFromList.append(distFromDigit)
                                    modalityList.append(modality)
                                    timePointList.append(h)
                                    trialList.append(k)
                                    runList.append(base)
                                    stimList.append(stimDigit)
                                    valList.append(val)
                                    distValList.append(j)
                                    layerList.append(layerNames[c])


data = pd.DataFrame({'subject': subList,
                     'volume': timePointList,
                     'modality': modalityList,
                     'data': valList,
                     'distFrom': distFromList,
                     'distVal': distValList,
                     'stim': stimList,
                     'run': runList,
                     'trialnr': trialList,
                     'layer': layerList})


# ===========================================================================
# Plotting

plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

for modality in ['BOLD']:
    tmp = data.loc[(data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]
    if modality == 'VASO':
        tmp = data.loc[(data['modality'] == modality) & (data['subject'] != 'sub-05') & (data['subject'] != 'sub-10') & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]

    g = sns.FacetGrid(tmp, col="distVal", hue='layer')
    g.map(sns.lineplot, "volume", "data")

    # for i in range(3):
    g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

    for j in range(7):
        # g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
        # g.axes[1, j].set_title("")  # Remove titles in lower row
        # g.axes[2, j].set_title("")  # Remove titles in lower row

        # g.axes[3, j].set_xticks(ticks, stimVolIds)

        g.axes[0, j].axvspan(14, 30, color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
        g.axes[0, j].axhline(y=0, linestyle="--", color='w')
        if j == 0:
            g.axes[0, j].set_title(f"Distance from peak: {j * 2} mm")  # Remove titles in lower row
        if j != 0:
            g.axes[0, j].set_title(f"{j*2} mm")  # Remove titles in lower row

    g.add_legend()
    g.tight_layout()
    g.savefig(f'results/group_centerSurround_summary_{modality}.png',
              bbox_inches="tight",
              dpi=400)


for sub in data['subject'].unique():
    for modality in ['BOLD']:
        tmp = data.loc[(data['subject'] == sub) & (data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]
        if modality == 'VASO':
            tmp = data.loc[(data['modality'] == modality) & (data['subject'] != 'sub-05') & (data['subject'] != 'sub-10') & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]

        g = sns.FacetGrid(tmp, col="distVal", hue='layer')
        g.map(sns.lineplot, "volume", "data")

        # for i in range(3):
        g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

        for j in range(7):
            # g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
            # g.axes[1, j].set_title("")  # Remove titles in lower row
            # g.axes[2, j].set_title("")  # Remove titles in lower row

            # g.axes[3, j].set_xticks(ticks, stimVolIds)

            g.axes[0, j].axvspan(14, 30, color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
            g.axes[0, j].axhline(y=0, linestyle="--", color='w')
            if j == 0:
                g.axes[0, j].set_title(f"Distance from peak: {j * 2} mm")  # Remove titles in lower row
            if j != 0:
                g.axes[0, j].set_title(f"{j*2} mm")  # Remove titles in lower row

        g.add_legend()
        g.tight_layout()
        g.savefig(f'results/{sub}_centerSurround_summary_{modality}.png',
                  bbox_inches="tight",
                  dpi=400)


# ===========================================================================
# Get layer-specific peak and trough ratio

distanceList = []
ratioList = []
layerList = []
minList = []
maxList = []

for modality in ['BOLD']:
    for layer in data['layer'].unique():
        for dist in range(2, 7):

            tmp = data.loc[(data['modality'] == modality) &
                           (data['distVal'] == dist) &
                           (data['distFrom'] == (data['stim'])) &
                           (data['layer'] == layer) &
                           (data['volume'] >= 15) &
                           (data['volume'] <= 20)]
            maxTp = 0
            maxVal = 0
            minTp = 0
            minVal = 100

            for tp in tmp['volume'].unique():

                tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                if tmpVal >= maxVal:
                    maxVal = tmpVal
                    maxTp = tp
                if tmpVal <= minVal:
                    minVal = tmpVal
                    minTp = tp

            print(maxVal)
            print(minVal)

            ratio = maxVal / minVal
            ratio = minVal / maxVal

            distanceList.append(dist)
            ratioList.append(ratio)
            layerList.append(layer)
            minList.append(minVal)
            maxList.append(maxVal)

ratioData = pd.DataFrame({'ratio': ratioList,
                          'distance': distanceList,
                          'layer': layerList,
                          'min': minList,
                          'max': maxList})


layerNames = ['Deep', 'Middle', 'Superficial']

# Plot ratio per distance
fig, ax = plt.subplots()
for i, value in enumerate(ratioData.layer.unique()):
    ax = sns.regplot(x="distance", y="ratio", ax=ax,
                     data=ratioData[ratioData.layer == value],
                     label=layerNames[i], x_jitter=.15, ci=None)


xticks = np.arange(2, 7)
xlabels = xticks * 2
ax.set_xticks(xticks, xlabels, fontsize=12)

ylabels = np.arange(-8, 1, 2)
yticks = np.arange(ylabels[0], 1, 2)
ax.set_yticks(yticks, ylabels, fontsize=12)

ax.set_xlabel(f'Distance from peak [mm]', fontsize=14)
ax.set_ylabel(f'Minimum / maximum deflection [a.u.]', fontsize=14)

plt.legend()
plt.tight_layout()
plt.savefig(f'results/group_inhibitionRatio.png', bbox_inches="tight")
plt.show()


# ===========================================================================
# Get layer-specific absolute peak and trough values

distanceList = []
layerList = []
valList = []
typeList = []

for modality in ['BOLD']:
    for layer in data['layer'].unique():
        for dist in range(2, 7):

            tmp = data.loc[(data['modality'] == modality) &
                           (data['distVal'] == dist) &
                           (data['distFrom'] == (data['stim'])) &
                           (data['layer'] == layer) &
                           (data['volume'] >= 15) &
                           (data['volume'] <= 20)]
            maxTp = 0
            maxVal = 0
            minTp = 0
            minVal = 100

            for tp in tmp['volume'].unique():

                tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                if tmpVal >= maxVal:
                    maxVal = tmpVal
                    maxTp = tp
                if tmpVal <= minVal:
                    minVal = tmpVal
                    minTp = tp

            distanceList.append(dist)
            typeList.append("Trough")
            valList.append(minVal)
            layerList.append(layer)

            distanceList.append(dist)
            typeList.append("Peak")
            valList.append(maxVal)
            layerList.append(layer)


peaksData = pd.DataFrame({'data': valList,
                          'distance': distanceList,
                          'layer': layerList,
                          'type': typeList
                          })


layerNames = ['Deep', 'Middle', 'Superficial']

g = sns.FacetGrid(peaksData, col="layer", hue='type')
g.map(sns.barplot, "distance", "data")


# Plot change of peak and trough independently
fig, ax = plt.subplots(1, 3)
for i, value in enumerate(ratioData.layer.unique()):
    tmp = peaksData.loc[peaksData['layer'] == value]
    sns.catplot(data=tmp, kind="bar", x="distance", y="val", hue="type", hue_order=['Peak', 'Trough'], ax=ax[i])

# xticks = np.arange(2, 7)
# xlabels = xticks * 2
# ax.set_xticks(xticks, xlabels, fontsize=12)
#
# ylabels = np.arange(-0.5, 0.6, 0.2).round(decimals=1)
# yticks = np.arange(-0.5, 0.6, 0.2).round(decimals=1)
# ax.set_yticks(yticks, ylabels, fontsize=12)
#
# ax.set_xlabel(f'Distance from peak [mm]', fontsize=14)
# ax.set_ylabel(f'Signal change [%]', fontsize=14)
#
# ax.axhline(y=0, linestyle="--", color='w')

plt.legend()
plt.tight_layout()
# plt.savefig(f'results/group_peakTrough_absolute.png', bbox_inches="tight")
plt.show()
