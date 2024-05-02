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
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

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

layerCompartments = [[1, 4], [5, 7], [8, 11]]

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

    for modality in ['VASO', 'BOLD']:

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

                for k, layerCompartment in enumerate(layerCompartments):

                    mask_4 = idxLayers >= layerCompartment[0]
                    mask_5 = idxLayers <= layerCompartment[1]
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

                    np.savetxt(f'results/{base}_distFrom{digit}_layer{k+1}.csv', timecourses, delimiter=',')


# ===========================================================================
# Extract distance dependent ERAs

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

    # for modality in ['BOLD']:
    for modality in ['VASO', 'BOLD']:
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

                for c in [1, 2, 3]:

                    # ====================================================
                    # Load previously extracted data
                    print('Load data')
                    data = np.loadtxt(f'results/{base}_distFrom{distFromDigit}_layer{c}.csv', delimiter=',')
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
                                    layerList.append(layerNames[c-1])


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

data.to_csv(f'results/distanceFromPeak.csv',
            sep=',',
            index=False)

# ===========================================================================
# Plotting

data = pd.read_csv(f'results/distanceFromPeak.csv', sep=',')
data.columns


plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

colors = {'#003f5c': ['#003f5c', '#8699aa', '#ffffff'],
          '#374c80': ['#374c80', '#9ba1be', '#ffffff'],
          '#7a5195': ['#7a5195', '#bca5c9', '#ffffff'],
          '#bc5090': ['#bc5090', '#e1a9c6', '#ffffff'],
          '#ef5675': ['#ef5675', '#ffafb7', '#ffffff'],
          '#ff764a': ['#ff764a', '#ffbca2', '#ffffff'],
          '#ffa600': ['#ffa600', '#ffd291', '#ffffff']
          }
hexes = ['#003f5c', '#374c80', '#7a5195', '#bc5090', '#ef5675', '#ff764a', '#ffa600']

for k, color in enumerate(hexes):

    for modality in ['BOLD']:
        tmp = data.loc[(data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]
        if modality == 'VASO':
            tmp = data.loc[(data['modality'] == modality) & (data['subject'] != 'sub-05') & (data['subject'] != 'sub-10') & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]

        g = sns.FacetGrid(tmp, col="distVal", hue='layer', palette=colors[color])
        g.map(sns.lineplot, "volume", "data")

        # for i in range(3):
        g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

        for j in range(7):
            # g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
            # g.axes[1, j].set_title("")  # Remove titles in lower row
            # g.axes[2, j].set_title("")  # Remove titles in lower row

            # g.axes[3, j].set_xticks(ticks, stimVolIds)
            g.axes[0, j].axvspan(13, 13 + (30/TR), color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
            g.axes[0, j].axhline(y=0, linestyle="--", color='w')

            g.axes[0, j].set_title(f"{j*2} - {(j+1)*2} mm", color=hexes[j])  # Remove titles in lower row

        g.add_legend()
        g.tight_layout()
        g.savefig(f'results/group_centerSurround_summary_{modality}_{k}.png',
                  bbox_inches="tight",
                  dpi=400)
        # g.close()

for sub in subs:
    for k, color in enumerate(hexes):
        for modality in ['BOLD', 'VASO']:
            tmp = data.loc[(data['subject'] == sub) & (data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]
            if modality == 'VASO':
                tmp = data.loc[(data['subject'] == sub) & (data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == (data['stim']))]

            g = sns.FacetGrid(tmp, col="distVal", hue='layer', palette=colors[color])
            g.map(sns.lineplot, "volume", "data")

            # for i in range(3):
            g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

            for j in range(7):
                # g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
                # g.axes[1, j].set_title("")  # Remove titles in lower row
                # g.axes[2, j].set_title("")  # Remove titles in lower row

                # g.axes[3, j].set_xticks(ticks, stimVolIds)

                g.axes[0, j].axvspan(13, 13 + (30 / TR), color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
                g.axes[0, j].axhline(y=0, linestyle="--", color='w')

                g.axes[0, j].set_title(f"{j*2} - {(j+1)*2} mm", color=hexes[j])  # Remove titles in lower row

            g.add_legend()
            g.tight_layout()
            g.savefig(f'results/{sub}_centerSurround_summary_{modality}_{k}.png',
                      bbox_inches="tight",
                      dpi=400)

# ==========================================================================================
# Check pattern of individual digits

for k, color in enumerate(hexes):
    for sub in ['sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17']:
        for modality in ['BOLD']:
            for digit in digits:
                tmp = data.loc[(data['subject'] == sub) & (data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == digit) & (data['stim'] == digit)]

                distVals = np.max(tmp['distVal'])

                g = sns.FacetGrid(tmp, col="distVal", hue='layer', palette=colors[color])
                g.map(sns.lineplot, "volume", "data")

                # for i in range(3):
                g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

                for j in range(distVals+1):
                    g.axes[0, j].axvspan(13, 13 + (30 / TR), color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
                    g.axes[0, j].axhline(y=0, linestyle="--", color='w')

                    g.axes[0, j].set_title(f"{j*2} - {(j+1)*2} mm", color=hexes[j])  # Remove titles in lower row

                g.add_legend()
                g.tight_layout()
                g.savefig(f'results/{sub}_centerSurround_summary_{modality}_{k}_{digit}.png',
                          bbox_inches="tight",
                          dpi=400)


for k, color in enumerate(hexes):
    for modality in ['VASO']:
        tmp = data.loc[(data['subject'] != 'sub-05') & (data['subject'] != 'sub-06') & (data['subject'] != 'sub-07') & (data['subject'] != 'sub-09') & (data['subject'] != 'sub-10') & (data['subject'] != 'sub-18') & (data['modality'] == modality) & (data['distVal'] <= 6) & (data['distFrom'] == data['stim'])]

        distVals = np.max(tmp['distVal'])

        g = sns.FacetGrid(tmp, col="distVal", hue='layer', palette=colors[color])
        g.map(sns.lineplot, "volume", "data")

        # for i in range(3):
        g.axes[0, 0].set_ylabel(f'Signal change [%]', fontsize=14)

        for j in range(distVals+1):
            g.axes[0, j].axvspan(13, 13 + (30/TR), color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
            g.axes[0, j].axhline(y=0, linestyle="--", color='w')

            g.axes[0, j].set_title(f"{j*2} - {(j+1)*2} mm", color=hexes[j])  # Remove titles in lower row

        g.add_legend()
        g.tight_layout()
        g.savefig(f'results/subGroup_centerSurround_summary_{modality}_{k}.png',
                  bbox_inches="tight",
                  dpi=400)


# Plot more distant and closest distances in one plot
TR = 1.9295

names = ['0-2 mm', '10-12 mm']

for modality in ['BOLD', 'VASO']:
    fig, ax = plt.subplots()

    for j, distbin in enumerate([0, 3]):
        if distbin == 0:
            tmp = data.loc[(data['modality'] == modality) & (data['distVal'] == distbin) & (data['distFrom'] == data['stim'])]
        if distbin != 0:
            tmp = data.loc[
                (data['modality'] == modality) & (data['distVal'] >= distbin) & (data['distVal'] < 6) & (data['distFrom'] == data['stim'])]

        sns.lineplot(data=tmp, y='data', x='volume', linewidth=2, color=colors[hexes[distbin]][1])

    ax.set_ylabel(f'{modality}\nSignal change [%]', fontsize=18)
    ax.axhline(y=0, linestyle="--", color='w')
    ax.axvspan(13, 13 + (30 / TR), color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')

    if modality == 'BOLD':
        plt.yticks(np.arange(-1, 2.5, 0.5), fontsize=14)

    if modality == 'VASO':
        plt.yticks(np.arange(-0.5, 1.1, 0.5), fontsize=14)

    plt.xticks(np.arange(0, np.max(tmp['volume'].unique())+1, 10), fontsize=14)
    ax.set_xlabel(f'Volume', fontsize=18)
    plt.tight_layout()
    plt.savefig(f'results/group_{modality}_0vs10mmdist.png', bbox_inches="tight")
    plt.show()



ticks = np.arange(4, 20, 2)
plt.xticks(ticks, fontsize=18)
# plt.xlim([0, 50])

plt.xlabel(f'# Censored volumes', fontsize=18)

plt.legend(title='Modality', loc='upper right', labels=['Notnulled', 'Nulled'], fontsize=18, title_fontsize=20)
plt.tight_layout()
plt.savefig(f'results/group_scrubbedVols.png', bbox_inches="tight")
plt.show()
