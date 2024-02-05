"""

Extracts statiscical values across cortical depth from digit ROIs.

"""

import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Set data path
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10']
subs = ['sub-05', 'sub-06']

# Set sessions to work on
sessions = ['ses-02']
ses = 'ses-01'
DIGITS = ['D2', 'D3', 'D4']

MODALITIES = ['BOLD', 'VASO']

subList = []
depthList = []
valList = []
roiList = []
stimList = []
modalityList = []
contrastList = []
roiModalityList = []

for sub in subs:
    print(f'Processing {sub}')

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    # make folder to dump statistocal maps
    statFolder = f'{funcDir}/statMaps'
    roiFolder = f'{funcDir}/rois'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    depthFile = f'{anatFolder}/seg_rim_polished_layers_equivol.nii.gz'
    depthNii = nb.load(depthFile)
    depthData = depthNii.get_fdata()
    layers = np.unique(depthData)[1:]

    for roiModality in ['VASO', 'BOLD']:
        roisData = nb.load(f'{roiFolder}/{sub}_{roiModality}_allRois.nii.gz').get_fdata()

        for i, roi in enumerate(DIGITS, start=1):
            print(f'Extracting from {roi} ROI')
            # roiIdx = roisData[roisData == i].astype('bool')
            roiIdx = roisData == i

            for digitStim in DIGITS:
                print(f'Extracting {digitStim} stimulation')

                for modality in MODALITIES:
                    print(f'{modality}')

                    for contrast in ['vsAll', 'vsRest']:

                        stimData = nb.load(f'{statFolder}/{sub}_{digitStim}{contrast}_{modality}_registered.nii').get_fdata()

                        # stimData = stimData[roiIdx.astype('bool')]
                        # depth = depthData[roiIdx.astype('bool')]

                        # for metric, val in zip(depth,stimData):
                        for layer in layers:

                            layerIdx = depthData == layer
                            tmp = roiIdx*layerIdx

                            val = np.mean(stimData[tmp])

                            subList.append(sub)
                            depthList.append(layer)
                            valList.append(val)
                            roiList.append(roi)
                            stimList.append(digitStim)
                            modalityList.append(modality)
                            contrastList.append(contrast)
                            roiModalityList.append(roiModality)


data = pd.DataFrame({'subject': subList,
                     'depth': depthList,
                     'value': valList,
                     'roi': roiList,
                     'stim': stimList,
                     'modality': modalityList,
                     'contrast': contrastList,
                     'roiModality': roiModalityList
                     })

data.to_csv(f'results/groupStimulationDepthData.csv', index=False)

# ==========================================
# Plotting
# ==========================================
plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

for contrast in ['vsAll', 'vsRest']:
    for sub in subs:
        for modality in ['BOLD', 'VASO']:
            for roi in DIGITS:
                tmp = data.loc[(data['subject'] == sub) & (data['modality'] == modality) & (data['roi'] == roi) & (data['contrast'] == contrast)]
                sns.lineplot(data=tmp, x='depth', y='value', hue='stim')
                plt.ylabel(f'z-score', fontsize=24)
                plt.title(f"{sub} {roi}-ROI {modality}", fontsize=24, pad=20)
                # plt.xlabel('WM                                CSF', fontsize=24)
                plt.xticks([])
                # yLimits = ax.get_ylim()
                # plt.ylim(0,yLimits[1])

                plt.yticks(fontsize=18)

                # ax.yaxis.set_major_locator(MaxNLocator(integer=True))

                plt.legend(loc='upper left')

                plt.savefig(f'results/{sub}_{roi}_{contrast}_zScoreProfile_{modality}.png', bbox_inches="tight")
                plt.show()

# ========================================
# Summary plot

for contrast in ['vsAll', 'vsRest']:
    tmp = data.loc[data['contrast'] == contrast]
    g = sns.FacetGrid(tmp, col="roi",  row="modality", hue='stim', sharey=False)
    g.map(sns.lineplot, "depth", "value")
    g.add_legend()
    plt.savefig(f'results/group_zScoreProfile_{contrast}_summary.png', bbox_inches="tight")
    plt.show()

# subs = ['sub-14']

for sub in subs:
    for roiModality in ['VASO', 'BOLD']:
        for contrast in ['vsAll', 'vsRest']:
            tmp = data.loc[(data['subject'] == sub) & (data['contrast'] == contrast) & (data['roiModality'] == roiModality)]
            g = sns.FacetGrid(tmp, col="roi",  row="modality", hue='stim', sharey=False)
            g.map(sns.lineplot, "depth", "value")
            g.add_legend()
            plt.savefig(f'results/{sub}_zScoreProfile_{contrast}_roi-{roiModality}_summary.png', bbox_inches="tight")
            plt.show()

