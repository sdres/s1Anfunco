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
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12']
# subs = ['sub-16', 'sub-17', 'sub-18']

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

    # for roiModality in ['VASO', 'BOLD']:
    for roiModality in ['BOLD']:
        roisData = nb.load(f'{roiFolder}/{sub}_{roiModality}_allRois.nii.gz').get_fdata()

        for i, roi in enumerate(DIGITS, start=1):
            print(f'Extracting from {roi} ROI')
            roiIdx = roisData == i

            for digitStim in DIGITS:
                print(f'Extracting {digitStim} stimulation')

                for modality in MODALITIES:
                    print(f'{modality}')

                    for contrast in ['vsAll', 'vsRest']:
                        statsFile = f'{statFolder}/{sub}_{digitStim}{contrast}_{modality}_registered.nii'
                        stimData = nb.load(statsFile).get_fdata()

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
data = pd.read_csv(f'results/groupStimulationDepthData.csv')

# ==========================================
# Plotting
# ==========================================
plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

# Define figzize
FS = (8, 5)
# define linewidth to 2
LW = 2
# Define fontsize size for x- and y-labels
labelSize = 24
# Define fontsize size for x- and y-ticks
tickLabelSize = 18
# Define fontsize legend text
legendTextSize = 18

# ========================================
# Summary plot

ticks = {'BOLD': range(-6, 12, 2),
         'VASO': range(-2, 4, 1)
         }

# for contrast in ['vsRest']:
for contrast in ['vsAll', 'vsRest']:
    for roiModality in ['VASO', 'BOLD']:
    # for roiModality in ['BOLD']:

        if roiModality == 'BOLD':
            tmp = data.loc[(data['contrast'] == contrast) & (data['roiModality'] == 'BOLD')]

        if roiModality == 'VASO':
            tmp = data.loc[(data['contrast'] == contrast)
                           & (data['roiModality'] == 'BOLD')
                           & (data['subject'] != 'sub-05')
                           & (data['subject'] != 'sub-10')]

        g = sns.FacetGrid(tmp, col="roi",  row="modality", hue='stim', sharey="row", palette=paletteDigits)
        g.map(sns.lineplot, "depth", "value")

        # Adapt y-axis and labels
        for i, modality in zip(range(2), ['BOLD', 'VASO']):
            g.axes[i, 0].set_ylabel(f"{modality}\nz-score", fontsize=14)
            g.axes[i, 0].set_yticks(ticks[modality], ticks[modality], fontsize=12)

        # Adapt x-axis and labels
        for j in range(3):
            g.axes[1, j].set_xticks([1, 11], ['WM', 'CSF'], fontsize=12)
            g.axes[1, j].set_xlabel('', fontsize=12)

        # Set titles and color spines
        for j, digit in zip(range(3), DIGITS):
            g.axes[0, j].set_title(f"{digit} ROI", fontsize=14, color=paletteDigits[digit])  # Set titles
            g.axes[1, j].set_title("")  # Remove titles in lower row
            # # Color spines
            if contrast == 'vsRest':
                for i in range(2):
                    g.axes[i, j].axhline(y=0, linestyle="--", color='w')
                #     g.axes[i, j].tick_params(color=paletteDigits[digit])
                #     for spine in g.axes[i, j].spines.values():
                #         spine.set_edgecolor(paletteDigits[digit])

        plt.tight_layout()
        plt.savefig(f'results/group_zScoreProfile_{contrast}_roi-{roiModality}_summary.png',
                    bbox_inches="tight",
                    dpi=400)
        plt.show()


# Single subs
for sub in subs:
    for contrast in ['vsAll', 'vsRest']:
        # for roiModality in ['VASO', 'BOLD']:
        for roiModality in ['BOLD']:
            tmp = data.loc[(data['subject'] == sub) & (data['contrast'] == contrast) & (data['roiModality'] == roiModality)]
            g = sns.FacetGrid(tmp, col="roi", row="modality", hue='stim', sharey="row", palette=paletteDigits)
            g.map(sns.lineplot, "depth", "value")

            # Adapt y-axis and labels
            for i, modality in zip(range(2), ['BOLD', 'VASO']):
                g.axes[i, 0].set_ylabel(f"{modality}\nz-score", fontsize=14)
                g.axes[i, 0].set_yticks(ticks[modality], ticks[modality], fontsize=12)

            # Adapt x-axis and labels
            for j in range(3):
                g.axes[1, j].set_xticks([1, 11], ['WM', 'CSF'], fontsize=12)
                g.axes[1, j].set_xlabel('', fontsize=12)

            # Set titles and color spines
            for j, digit in zip(range(3), DIGITS):
                g.axes[0, j].set_title(f"{digit} ROI", fontsize=14, color=paletteDigits[digit])  # Set titles
                g.axes[1, j].set_title("")  # Remove titles in lower row

                if contrast == 'vsRest':
                    for i in range(2):
                        g.axes[i, j].axhline(y=0, linestyle="--", color='w')

                # # Color spines
                # for i in range(2):
                #     g.axes[i, j].tick_params(color=paletteDigits[digit])
                #     for spine in g.axes[i, j].spines.values():
                #         spine.set_edgecolor(paletteDigits[digit])

            plt.tight_layout()
            plt.savefig(f'results/{sub}_zScoreProfile_{contrast}_roi-{roiModality}_summary.png',
                        bbox_inches="tight",
                        dpi=400)

            # plt.show()
            plt.close()

