"""Extract cortical thinkness from perimeter"""

import nibabel as nb
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
DIGITS = ['D2', 'D3', 'D4']

ROOT = f'/Users/sebastiandresbach/data/s1Anfunco/Nifti'

voxelList = []
subList = []
thicknessList = []
roiModalityList = []
roiList = []

for sub in subs:
    rimFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'
    funcDir = f'{ROOT}/derivatives/{sub}/func'
    roiFolder = f'{funcDir}/rois'

    thickFile = f'{rimFolder}/seg_rim_polished_thickness.nii.gz'
    thickNii = nb.load(thickFile)
    thickData = thickNii.get_fdata()

    for roiModality in ['BOLD']:
        roisData = nb.load(f'{roiFolder}/{sub}_{roiModality}_allRois.nii.gz').get_fdata()
        for i, roi in enumerate(DIGITS, start=1):
            roiIdx = roisData == i

            thickness = thickData[roiIdx]

            for j, thick in enumerate(thickness):
                subList.append(sub)
                voxelList.append(j)
                thicknessList.append(thick)
                roiModalityList.append(roiModality)
                roiList.append(roi)

data = pd.DataFrame({'subject': subList,
                     'idx': voxelList,
                     'roi': roiList,
                     'thickness': thicknessList,
                     'roiModality': roiModalityList})


plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}

mean = np.mean(data['thickness'])

fig, ax = plt.subplots(figsize=(8, 5))
sns.boxplot(data, ax=ax, x='subject', y='thickness', showfliers=False)
ax.axhline(y=mean, linestyle="--", color='w')
ax.set_title('Gray matter ROI thickness', fontsize=16)
ax.set_xticks(range(len(subs)), subs, fontsize=12)
ax.set_xlabel('', fontsize=12)
ax.set_yticks(range(1, 5), range(1, 5), fontsize=12)
ax.set_ylabel('Cortical thickness [mm]', fontsize=14)
plt.tight_layout()
plt.savefig(f'results/group_cortThick.png')
plt.show()
