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
# subs = ['sub-05']
# subs = ['sub-17', 'sub-18']
skipVols = 2

TR = 1.9295

# ===========================================================================
# Limit timecourse to peri

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regStimDir = f'{funcDir}/registeredStim'

    # for modality in ['BOLD']:
    for modality in ['VASO', 'BOLD']:
        print(f'Processing {modality}')

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)

            # ================================================
            # Load functional data and combine into timeseries

            volumes = sorted(glob.glob(f'{regStimDir}/peri/{base}*vol*_registered.nii.gz'))

            vol1 = nb.load(volumes[0]).get_fdata()
            idx = vol1 > 0

            shape = vol1[idx].shape
            new = np.zeros((shape[0], len(volumes)), dtype="float16")

            for i, vol in enumerate(volumes):
                if i % 10 == 0:
                    print(f'Adding volume {i}')

                nii = nb.load(vol)
                tmp = np.asarray(nii.dataobj)[idx]

                new[:, i] = tmp

            np.savetxt(f'{regStimDir}/{base}_maskedPeri.csv', new[:, skipVols:], delimiter=',')


# ===========================================================================
# Clean timecourse

for sub in subs:
    print(f'Processing {sub}')
    funcDir = f'{ROOT}/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Make dir to temporarily dump registered volumes
    regStimDir = f'{funcDir}/registeredStim'

    # for modality in ['VASO']:
    for modality in ['VASO', 'BOLD']:
        print(f'Processing {modality}')

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)

            # ================================================
            # Load motion timecourse
            motionDataFile = glob.glob(f'{ROOT}/{sub}/func/motionParameters/{base[:-5]}/{base}_nuisance.txt')[0]

            motionData = pd.read_csv(motionDataFile, delimiter='\t', header=None)

            # ================================================
            # Merge regressors

            # motionData['csf'] = csfData

            regressors = motionData

            # ================================================
            # Load functional data

            new = np.loadtxt(f'{regStimDir}/{base}_maskedPeri.csv', delimiter=',')

            minVal = np.min(new)
            maxVal = np.max(new)
            meanVal = np.mean(new)
            print(f'Before cleaning')
            print(f'Mean: {meanVal}')
            print(f'Range: {minVal} - {maxVal}')
            print('')

            new = new.transpose()

            regressorsNew = regressors.iloc[skipVols:]

            cleaned = signal.clean(new,
                                   detrend=True,
                                   standardize='psc',
                                   # standardize=False,
                                   confounds=regressorsNew,
                                   high_pass=0.01,
                                   t_r=TR
                                   )

            np.savetxt(f'{regStimDir}/{base}_maskedPeri_clean_psc.csv', cleaned, delimiter=',')

#             minVal = np.min(cleaned)
#             maxVal = np.max(cleaned)
#             meanVal = np.mean(cleaned)
#             print(f'After cleaning')
#             print(f'Mean: {meanVal}')
#             print(f'Range: {minVal} - {maxVal}')
#             print('')


# ===========================================================================
# Check range
#
# for sub in subs:
#     print(f'Processing {sub}')
#     funcDir = f'{ROOT}/{sub}/func'
#     anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
#
#     # Make dir to temporarily dump registered volumes
#     regStimDir = f'{funcDir}/registeredStim'
#
#     for modality in ['BOLD']:
#     # for modality in ['VASO', 'BOLD']:
#         print(f'Processing {modality}')
#
#         # Find subject-specific stimulation run(s)
#         stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))
#
#         for stimRun in stimRuns:
#             base = os.path.basename(stimRun).rsplit('.', 2)[0]
#             print(base)
#
#             # ================================================
#             # Load functional data
#
#             new = np.loadtxt(f'{regStimDir}/{base}_maskedPeri.csv', delimiter=',')
#
#             minVal = np.min(new)
#             maxVal = np.max(new)
#             meanVal = np.mean(new)
#             print(f'Before cleaning')
#             print(f'Mean: {meanVal}')
#             print(f'Range: {minVal} - {maxVal}')
#             print('')
#
#             cleaned = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean.csv', delimiter=',')
#
#             minVal = np.min(cleaned)
#             maxVal = np.max(cleaned)
#             meanVal = np.mean(cleaned)
#             print(f'After cleaning')
#             print(f'Mean: {meanVal}')
#             print(f'Range: {minVal} - {maxVal}')
#             print('')
# #             np.argwhere(cleaned < 0)
# #
# test = cleaned.transpose()
# # test.shape
# #
# x = np.arange(test.shape[1])
# plt.plot(x, test[108099, :])
# plt.plot(x, new.transpose()[108099, :])

# ===========================================================================
# Calculate PSC
#
# digits = ['D2', 'D3', 'D4']
# # subs = ['sub-05']
# for sub in subs:
#     print(f'Processing {sub}')
#
#     # ====================================================
#     # Set folders
#     funcDir = f'{ROOT}/{sub}/func'
#     anatFolder = f'{ROOT}/{sub}/anat/upsampled'
#     regStimDir = f'{funcDir}/registeredStim'
#
#     for modality in ['BOLD']:
#     # for modality in ['VASO', 'BOLD']:
#         print(f'Extracting from {modality}')
#
#         # Find subject-specific stimulation run(s)
#         stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))
#
#         for stimRun in stimRuns:
#             base = os.path.basename(stimRun).rsplit('.', 2)[0]
#             print(base)
#
#             # ====================================================
#             # Load previously extracted data
#             data = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean.csv', delimiter=',')
#
#             data = data.transpose()
#
#             nrVols = data.shape[1]
#
#             # ====================================================
#             # Get baseline
#
#             # Find non stim-periods
#             # Load stim info
#             print('Load stimulation times')
#             stimInfo = pd.DataFrame()
#             for digit in digits:
#                 tmpFile = f'{ROOT}/designFiles/stimulationPatterns/{sub}/{base[:-5]}_{digit}.txt'
#                 tmp = pd.read_csv(tmpFile, names=['start', 'dur', 'mod'], sep=' ')
#
#                 stimInfo = pd.concat([stimInfo, tmp])
#
#             stimInfo = stimInfo.sort_values(by=['start'])
#             stimInfo = stimInfo.reset_index(drop=True)
#             # print(stimInfo)
#
#             # Get baseline volumes
#             print('Calculate baseline')
#             # Check for initial baseline
#             if stimInfo['start'].to_numpy()[0] <= 30:
#                 print(f'{sub} has no initial baseline')
#                 # Get end of last stimulation period
#                 lastStimEnd = (stimInfo['start'].to_numpy()[-1] + 60) / TR
#                 baselineVols = np.arange(lastStimEnd, nrVols).astype('int')
#
#             else:
#                 # Take first and final rest period
#                 baselineVols = np.arange(int(30 / TR))
#                 period2 = np.arange(nrVols - int(30 / TR), nrVols)
#
#                 baselineVols = np.concatenate([baselineVols, period2]).astype('int')
#
#             # for i, onset in enumerate(stimInfo['start']):
#             #     # include all volumes until first onset
#             #     if i == 0:
#             #         baselineVols = np.arange(int(onset / TR))
#             #
#             #     # For all other rest periods include the last 10 seconds
#             #     elif i != 0:
#             #         first = int((onset - 10) / TR)
#             #         last = int((onset / TR) + TR)
#             #         tmp = np.arange(first, last)
#             #         baselineVols = np.concatenate([baselineVols, tmp])
#             #
#             # # Include last rest period completely
#             # endLastStim = (stimInfo['start'].to_numpy()[-1] + stimInfo['dur'].to_numpy()[-1]) / TR
#             # tmp = np.arange(int(endLastStim), nrVols)
#             # baselineVols = np.concatenate([baselineVols, tmp])
#             # test = np.take(data, baselineVols, axis=1)
#             # test.shape
#             # baselineVols.shape
#             #
#             # x = np.arange(data.shape[1])
#             # plt.plot(x, data[100000])
#
#             # Calculate baseline
#             bl = np.mean(np.take(data, baselineVols, axis=1), axis=1)
#
#             # Calculate PSC
#             print('Calculate PSC')
#             psc = np.zeros(data.shape)
#
#             for j in range(data.shape[-1]):
#
#                 tmp = (np.divide(data[:, j], bl) - 1) * 100
#                 psc[:, j] = tmp
#
#             # minVal = np.min(psc)
#             # maxVal = np.max(psc)
#             # meanVal = np.mean(psc)
#             # print(f'After PSC calculation')
#             # print(f'Mean: {meanVal}')
#             # print(f'Range: {minVal} - {maxVal}')
#             # print('')
#             print('Save psc data')
#
#             # Calculate std dev
#             stdDevs = np.std(psc, axis=1)
#             # Remove stdDev higher than 10
#             counter = 0
#             for i, val in enumerate(stdDevs):
#                 if val >= 10:
#                     psc[i, :] = np.nan
#                     counter += 1
#             print(f'Removed data from {counter} voxels')
#
#             np.savetxt(f'{regStimDir}/{base}_maskedPeri_clean_psc.csv', psc, delimiter=',')


# ===========================================================================
# Extract timecourse

digits = ['D2', 'D3', 'D4']

for sub in ['sub-09']:
    print(f'Processing {sub}')

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

    for modality in ['VASO']:
    # for modality in ['VASO', 'BOLD']:
        print(f'Extracting from {modality}')

        # Find subject-specific stimulation run(s)
        stimRuns = sorted(glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}.nii.gz'))

        for stimRun in stimRuns:
            base = os.path.basename(stimRun).rsplit('.', 2)[0]
            print(base)

            # ====================================================
            # Load previously extracted data
            # data = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean.csv', delimiter=',')
            data = np.loadtxt(f'{regStimDir}/{base}_maskedPeri_clean_psc.csv', delimiter=',')

            for vox in range(data.shape[1]):
                tmp = data[:, vox]
                if not (tmp >= 5).any():
                    print(f'Threshold not exceeded')


            data.shape
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

            np.savetxt(f'results/{base}_layerTimecourses_clean_psc.csv', timecourses, delimiter=',')


# ===========================================================================
# Generate ERAs

digits = ['D2', 'D3', 'D4']

layerCompartments = {0: [1, 4],
                     1: [4, 7],
                     2: [7, 10]
                     }

layerNames = ['Deep', 'Middle', 'Superficial']

timePointList = []
modalityList = []
valList = []
stimList = []
subList = []
depthList = []
roiList = []
runList = []
trialList = []

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

            # ====================================================
            # Load previously extracted data
            print('Load data')
            data = np.loadtxt(f'results/{base}_layerTimecourses_clean_psc.csv', delimiter=',')
            data = data.transpose()

            offset = 0
            # Loop over ROIs
            for j, roiDigit in enumerate(digits):

                for c in layerCompartments:

                    firstLay = layerCompartments[c][0]
                    lastLay = layerCompartments[c][1]

                    # print(f"{roiDigit} layer compartment {c+1}")
                    # print(f'extracting from dimension {firstLay + offset} until {lastLay + offset}')

                    # Load stim info
                    for i, stimDigit in enumerate(digits):

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

                            tmpData = data[firstLay + offset: lastLay + offset, start: end]
                            tmpData = np.mean(tmpData, axis=0)

                            if modality == 'VASO':
                                tmpData = tmpData * -1

                            for h, val in enumerate(tmpData):

                                subList.append(sub)
                                depthList.append(layerNames[c])
                                roiList.append(roiDigit)
                                modalityList.append(modality)
                                timePointList.append(h)
                                trialList.append(k)
                                runList.append(base)
                                stimList.append(stimDigit)
                                valList.append(val)

                offset = offset + 11

data = pd.DataFrame({'subject': subList,
                     'volume': timePointList,
                     'modality': modalityList,
                     'data': valList,
                     'layer': depthList,
                     'roi': roiList,
                     'stim': stimList,
                     'run': runList,
                     'trialnr': trialList})


data.to_csv(f'results/ERAs.csv',
            sep=',',
            index=False)

plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

layerNames_rev = layerNames[::-1]

for sub in subs:
    for modality in ['BOLD', 'VASO']:
        tmp = data.loc[(data['modality'] == modality) & (data['subject'] == sub)]

        stimVolIds = np.arange(16)

        stimVolIds = np.concatenate([[-2, -1], stimVolIds])
        ticks = np.arange(len(stimVolIds))

        g = sns.FacetGrid(tmp, col="roi", row="layer", hue='stim', row_order=['Superficial', 'Middle', 'Deep'], palette=paletteDigits)
        g.map(sns.lineplot, "volume", "data")
        # For all axes
        for i in range(3):
            g.axes[i, 0].set_ylabel(f'{layerNames_rev[i]} layer\nSignal change [%]', fontsize=14)

            for j in range(3):
                g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
                g.axes[1, j].set_title("")  # Remove titles in lower row
                g.axes[2, j].set_title("")  # Remove titles in lower row

                # g.axes[3, j].set_xticks(ticks, stimVolIds)

                g.axes[i, j].axvspan(14, 30, color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
                g.axes[i, j].axhline(y=0, linestyle="--", color='w')

                # g.axes[i, j].set_ylim(-5, 10)

        g.add_legend()
        g.tight_layout()
        g.savefig(f'results/{sub}_ERA_{modality}_summary.png',
                    bbox_inches="tight",
                    dpi=400)
        # plt.show()


for modality in ['BOLD', 'VASO']:
# for modality in ['BOLD']:
    tmp = data.loc[(data['modality'] == modality) & (data['subject'] == sub)]
    tmp = data.loc[(data['modality'] == modality)]

    stimVolIds = np.arange(16)

    stimVolIds = np.concatenate([[-2, -1], stimVolIds])
    ticks = np.arange(len(stimVolIds))

    g = sns.FacetGrid(tmp, col="roi", row="layer", hue='stim', row_order=['Superficial', 'Middle', 'Deep'], palette=paletteDigits)
    g.map(sns.lineplot, "volume", "data")
    # For all axes
    for i in range(3):
        g.axes[i, 0].set_ylabel(f'{layerNames_rev[i]} layer\nSignal change [%]', fontsize=14)

        for j in range(3):
            g.axes[0, j].set_title(f"{digits[j]} ROI", fontsize=14, color=paletteDigits[digits[j]])  # Set titles
            g.axes[1, j].set_title("")  # Remove titles in lower row
            g.axes[2, j].set_title("")  # Remove titles in lower row

            # g.axes[3, j].set_xticks(ticks, stimVolIds)

            g.axes[i, j].axvspan(14, 30, color='#e5e5e5', alpha=0.2, lw=0, label='stimulation')
            g.axes[i, j].axhline(y=0, linestyle="--", color='w')


            # g.axes[i, j].set_ylim(-5, 10)

    g.add_legend()
    g.tight_layout()
    g.savefig(f'results/group_ERA_{modality}_summary.png',
                bbox_inches="tight",
                dpi=400)
    # plt.show()
