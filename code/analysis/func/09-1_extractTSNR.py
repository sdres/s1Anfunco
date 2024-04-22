"""

Extracts QA values from perimeter.

"""

import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob
import os

# Set data path
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10']
# subs = ['sub-06']

# Set sessions to work on
sessions = ['ses-02']
ses = 'ses-01'
DIGITS = ['D2', 'D3', 'D4']

MODALITIES = ['BOLD', 'VASO']

subList = []
valList = []
modalityList = []
voxList = []
metricList = []

for sub in subs:
    print(f'Processing {sub}')

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    perimeterFile = f'{anatFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    perimeterNii = nb.load(perimeterFile)
    perimeterData = perimeterNii.get_fdata()
    perimeterIdx = perimeterData == 1

    for modality in ['BOLD', 'VASO']:
        for metric in ['tSNR', 'mean']:
            statFiles = glob.glob(f'{funcDir}/{sub}_ses-0*_task-stim_run-00*_{modality}_{metric}_registered.nii')
            for statFile in statFiles:
                base = os.path.basename(statFile).rsplit('.', 2)[0][:-4]
                statNii = nb.load(statFile)
                statData = statNii.get_fdata()

                data = statData[perimeterIdx]

                for i, val in enumerate(data):

                    subList.append(sub)
                    voxList.append(i)
                    valList.append(val)
                    modalityList.append(modality)
                    metricList.append(metric)

data = pd.DataFrame({'subject': subList,
                     'voxel': voxList,
                     'value': valList,
                     'modality': modalityList,
                     'metric': metricList
                     })

data.to_csv(f'results/groupTsnrData.csv', index=False)

# ==========================================
# Plotting
# ==========================================
plt.style.use('dark_background')

data = pd.read_csv(f'results/groupTsnrData.csv')
data['subject'].unique()
palette = {
    'BOLD': 'tab:orange',
    'VASO': 'tab:blue'}

for metric in ['tSNR']:
    fig, ax = plt.subplots()

    tmp = data.loc[data['metric'] == metric]

    sns.kdeplot(data=tmp, x='value', hue='modality', linewidth=2, palette=palette)

    plt.title(f'Group (n=11)', fontsize=24)    # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylabel('Voxel count', fontsize=24)
    plt.yticks([])

    ticks = np.arange(0, 51, 10)
    plt.xticks(ticks, fontsize=18)
    plt.xlim([0, 50])

    plt.xlabel(f'{metric}', fontsize=24)

    #legend hack
    old_legend = ax.legend_
    handles = old_legend.legend_handles
    labels = ['BOLD', 'VASO']
    title = old_legend.get_title().get_text()
    ax.legend(handles, labels, loc='upper right', title='', fontsize=20)
    plt.tight_layout()
    plt.savefig(f'results/QA/group_{metric}.png', bbox_inches="tight")
    plt.close()


# Plot single participant data
for sub in subs:
    for metric in ['tSNR']:
        fig, ax = plt.subplots()

        tmp = data.loc[(data['metric'] == metric) & (data['subject'] == sub)]

        sns.kdeplot(data=tmp, x='value', hue='modality', linewidth=2, palette=palette)

        plt.title(f'{sub}', fontsize=24)
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.ylabel('Voxel count', fontsize=24)
        plt.yticks([])

        ticks = np.arange(0, 51, 10)
        plt.xticks(ticks, fontsize=18)
        plt.xlim([0, 50])

        plt.xlabel(f'{metric}', fontsize=24)
        ax.get_legend().remove()
        #legend hack
        # old_legend = ax.legend_
        # handles = old_legend.legend_handles
        # labels = ['BOLD', 'VASO']
        # title = old_legend.get_title().get_text()
        # ax.legend(handles, labels, loc='upper right', title='', fontsize=20)
        plt.tight_layout()
        plt.savefig(f'results/QA/{sub}_{metric}.png', bbox_inches="tight")

        plt.close()

# ==============================
# Check for VASO voxels > 1

meanVals = data.loc[(data['metric'] == 'mean') & (data['modality'] == 'VASO')]
inflow = meanVals.loc[meanVals['value'] > 100]

PALETTE = {'BOLD': 'tab:orange', 'VASO': 'tab:blue'}

tmp = meanVals.loc[meanVals['value'] > 100]

fig, ax = plt.subplots(1, 1, figsize=(7.5, 5))
sns.histplot(data=tmp, x='value', hue='subject', multiple='dodge', linewidth=1, bins=10)
plt.ylabel('# Voxels', fontsize=24)
plt.yticks(np.arange(0, 301, 50), fontsize=18)

ticks = np.arange(100, 131, 5)
plt.xticks(ticks, fontsize=18)
# plt.xlim([0, 50])

plt.xlabel(f'VASO voxel mean', fontsize=18)

plt.tight_layout()
plt.savefig(f'results/group_VASO_mean_greater1.png', bbox_inches="tight")
plt.show()


