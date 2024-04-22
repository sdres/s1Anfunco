
import nibabel as nb
import numpy as np
import subprocess
import os
import glob
from nilearn import signal
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# subs = ['sub-16']
# subs = ['sub-14', 'sub-15', 'sub-17', 'sub-18']

# subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12']

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'
digits = ['D2', 'D3', 'D4']


# Excluce ROIs

layerCompartments = {0: [1, 4],
                     1: [4, 7],
                     2: [7, 10]
                     }

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

# subs = ['sub-15']
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

    # ====================================================
    # Load digit ROIs
    roiFile = f'{funcDir}/rois/{sub}_BOLD_allRois.nii.gz'
    roiNii = nb.load(roiFile)  # Load nifti
    roiData = roiNii.get_fdata()  # Load data as array
    idxRois = roiData[idx]

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

                    # Excluding ROIs
                    mask7 = np.where(idxRois >= 1, False, mask_6)

                    timecourses = np.zeros((nrVols, nrRegions))

                    for vol in range(nrVols):
                        tmp = data[:, vol]

                        offset = 0

                        for i, distance in enumerate(distances[1:]):
                            # Mask distances
                            mask_1 = distIdx >= distances[i]
                            mask_2 = distIdx <= distance
                            mask_3 = np.logical_and(mask_1, mask_2)

                            mask = np.multiply(mask_3, mask7)

                            val = np.mean(tmp[mask])

                            timecourses[vol, i] = val

                    np.savetxt(f'results/{base}_distFrom{digit}_layer{c+1}_excludeRois.csv', timecourses, delimiter=',')



# ===========================================================================
# Extract distance dependent timecourses

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
                    data = np.loadtxt(f'results/{base}_distFrom{distFromDigit}_layer{c+1}_excludeRois.csv', delimiter=',')
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

            g.axes[0, j].axvspan(14, 30, color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
            g.axes[0, j].axhline(y=0, linestyle="--", color='w')

            g.axes[0, j].set_title(f"{j*2} - {(j+1)*2} mm", color=hexes[j])  # Remove titles in lower row

        g.add_legend()
        g.tight_layout()
        g.savefig(f'results/group_centerSurround_summary_{modality}_{k}_accountRois.png',
                  bbox_inches="tight",
                  dpi=400)
