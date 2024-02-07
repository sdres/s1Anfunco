"""

Extracts QA values from perimeter.

"""

import numpy as np
import nibabel as nb
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob

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
        for metric in ['tSNR']:
            statFile = glob.glob(f'{funcDir}/{sub}_ses-0*_task-restBA3b_run-001_{modality}_{metric}_registered.nii')[0]
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


palette = {
    'BOLD': 'tab:orange',
    'VASO': 'tab:blue'}

for metric in data['metric'].unique():
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
    plt.show()


# Plot single participant data
for sub in subs:
    for metric in data['metric'].unique():
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

        plt.show()
