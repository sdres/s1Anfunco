"""Extracting ERAs from ROIs"""

import nibabel as nb
import numpy as np
import os
import glob
from nilearn import signal
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'

# Set subjects to work on
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

skipVols = 2
digits = ['D2', 'D3', 'D4']

TR = 1.9295


# ===========================================================================
# Determine runs with Pattern A

layerNames = ['Deep', 'Middle', 'Superficial']

runList = []

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for sub in subs:
    print(f'Processing {sub}')

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'
    regStimDir = f'{funcDir}/registeredStim'

    # Find subject-specific stimulation run(s)
    stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_BOLD.nii.gz'))

    for stimRun in stimRuns:
        base = os.path.basename(stimRun).rsplit('.', 2)[0]
        print(base)
        runList.append(base)
        stimInfo = pd.DataFrame()
        fingersList = []
        for digit in digits:
            tmpFile = f'{ROOT}/designFiles/stimulationPatterns/{sub}/{base[:-5]}_{digit}.txt'
            tmp = pd.read_csv(tmpFile, names=['start', 'dur', 'mod'], sep=' ')

            for i in range(len(tmp)):
                fingersList.append(int(digit[-1]))

            stimInfo = pd.concat([stimInfo, tmp])

        stimInfo['finger'] = fingersList
        stimInfo = stimInfo.sort_values(by=['start'])
        stimInfo = stimInfo.reset_index(drop=True)
        # print(stimInfo)
        if sub == 'sub-05':
            orders = np.array(stimInfo['finger'])
        if sub != 'sub-05':
            orders = np.vstack([orders, np.array(stimInfo['finger'])])

runNames = np.array(runList)

uniqueOrders = np.unique(orders, axis=0)
runIndices = {}

# find indices of sets of orders
for i, uniqueOrder in enumerate(uniqueOrders):
    count = 0
    runIndices[i] = []
    for j, order in enumerate(orders):
        if (uniqueOrder == order).all():
            count += 1
            runIndices[i].append(j)
    print(f'Found {count} occurences of {uniqueOrder}')


runs = np.take(runNames, runIndices[0])


# ===========================================================================
# Clean timecourse

for run in runs:
    print(f'Processing {run}')
    sub = run[:6]

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'
    regStimDir = f'{funcDir}/registeredStim'

    for modality in ['BOLD']:
        print(f'Extracting from {modality}')

        # Find subject-specific stimulation run(s)
        stimRun = sorted(glob.glob(f'{funcDir}/{run}.nii.gz'))[0]
        base = os.path.basename(stimRun).rsplit('.', 2)[0]
        print(base)

        # ================================================
        # Load functional data

        new = np.loadtxt(f'{regStimDir}/{base}_maskedPeri.csv', delimiter=',')

        new = new.transpose()

        psc = signal.clean(new,
                           detrend=True,
                           standardize='psc',
                           confounds=None,
                           high_pass=None,
                           t_r=TR
                           )

        np.savetxt(f'{regStimDir}/{base}_maskedPeri_psc_detrend_highPass.csv', psc, delimiter=',')


# ===========================================================================
# Extract timecourse

for run in runs:
    print(f'Processing {run}')
    sub = run[:6]

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'
    regStimDir = f'{funcDir}/registeredStim'

    # ====================================================
    # get shape for index
    volumes = sorted(glob.glob(f'{regStimDir}/peri/*BOLD*vol*_registered.nii.gz'))
    vol1 = nb.load(volumes[0]).get_fdata()
    idx = vol1 > 0
    idx.shape

    # ====================================================
    # Get depth information
    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    idxLayers = depthData[idx]
    layers = np.unique(idxLayers)
    idxLayers.shape

    # ====================================================
    # Load digit ROIs
    roiFile = f'{funcDir}/rois/{sub}_BOLD_allRois.nii.gz'
    roiNii = nb.load(roiFile)  # Load nifti
    roiData = roiNii.get_fdata()  # Load data as array
    idxRois = roiData[idx]
    idxRois.shape

    for modality in ['BOLD']:
    # for modality in ['VASO', 'BOLD']:
        print(f'Extracting from {modality}')

        stimRun = sorted(glob.glob(f'{funcDir}/{run}.nii.gz'))[0]
        base = os.path.basename(stimRun).rsplit('.', 2)[0]

        # ====================================================
        # Load previously extracted data
        # data = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean.csv', delimiter=',')
        data = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_psc_detrend_highPass.csv', delimiter=',')

        # data = data.transpose()
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

                    val = np.nanmean(volData[tmp])

                    timecourses[i, (offset+k)] = val

                offset = offset + len(layers)

        np.savetxt(f'results/{base}_layerTimecourses_psc_detrend_highPass.csv', timecourses, delimiter=',')

# ===========================================================================
# Generate run averages

layerCompartments = {0: [1, 4],
                     1: [4, 7],
                     2: [7, 10]
                     }

cols = []
offset = 0
for j in range(3):
    for c in layerCompartments:
        tmp = np.arange(layerCompartments[c][0], layerCompartments[c][1])
        for i in tmp:
            cols.append(i + offset)
    offset += 11

layerNames = ['Deep', 'Middle', 'Superficial']

timePointList = []
modalityList = []
valList = []
subList = []
depthList = []
roiList = []
runList = []

for run in runs:
    print(f'Processing {run}')
    sub = run[:6]

    # ====================================================
    # Set folders
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/{sub}/anat/upsampled'
    regStimDir = f'{funcDir}/registeredStim'

    for modality in ['BOLD']:
    # for modality in ['VASO', 'BOLD']:
        print(f'Extracting from {modality}')

        # Find subject-specific stimulation run(s)
        stimRun = sorted(glob.glob(f'{funcDir}/{run}.nii.gz'))

        # ====================================================
        # Load previously extracted data
        print('Load data')
        data = np.loadtxt(f'results/{run}_layerTimecourses_psc_detrend_highPass.csv', delimiter=',')
        data = data.transpose()

        offset = 0
        # Loop over ROIs
        for j, roiDigit in enumerate(digits):

            for c in layerCompartments:

                firstLay = layerCompartments[c][0]
                lastLay = layerCompartments[c][1]

                tmpData = data[firstLay + offset: lastLay + offset, :]
                tmpData = np.mean(tmpData, axis=0)

                if modality == 'VASO':
                    tmpData = tmpData * -1

                for h, val in enumerate(tmpData):

                    subList.append(sub)
                    depthList.append(layerNames[c])
                    roiList.append(roiDigit)
                    modalityList.append(modality)
                    timePointList.append(h)
                    runList.append(base)
                    valList.append(val)

            offset = offset + 11

data = pd.DataFrame({'subject': subList,
                     'volume': timePointList,
                     'modality': modalityList,
                     'data': valList,
                     'layer': depthList,
                     'roi': roiList,
                     'run': runList})


plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

stimInfo = pd.DataFrame()
fingersList = []
for digit in digits:
    tmpFile = f'{ROOT}/designFiles/stimulationPatterns/{sub}/{run[:-5]}_{digit}.txt'
    tmp = pd.read_csv(tmpFile, names=['start', 'dur', 'mod'], sep=' ')
    for i in range(len(tmp)):
        fingersList.append(digit)

    stimInfo = pd.concat([stimInfo, tmp])

stimInfo['finger'] = fingersList
stimInfo = stimInfo.sort_values(by=['start'])
stimInfo = stimInfo.reset_index(drop=True)


# for modality in ['BOLD', 'VASO']:
for modality in ['BOLD']:

    # tmp = data.loc[(data['modality'] == modality) & (data['subject'] == sub)]
    tmp = data.loc[(data['modality'] == modality)]

    # stimVolIds = np.arange(16)

    # stimVolIds = np.concatenate([[-2, -1], stimVolIds])
    # ticks = np.arange(len(stimVolIds))

    g = sns.FacetGrid(tmp, row="roi", height=3, aspect=4, hue='roi', palette=paletteDigits)
    g.map(sns.lineplot, "volume", "data")

    # For all axes
    for i in range(3):
        g.axes[i, 0].set_ylabel(f'{digits[i]} ROI\nsignal change [%]', fontsize=14)
        g.axes[i, 0].set_title("")  # Remove titles in lower row
        g.axes[i, 0].axhline(y=0, linestyle="--", color='w')

        for j, row in stimInfo.iterrows():
            color = paletteDigits[row['finger']]
            start = row['start'] / TR
            end = (row['start'] + row['dur']) / TR + 1
            g.axes[i, 0].axvspan(start, end, color=color, alpha=0.2, lw=0, label='stimulation')

    g.tight_layout()
    g.savefig(f'results/group_entireTS_{modality}_summary_detrend_highPass.png',
                bbox_inches="tight",
                dpi=400)
