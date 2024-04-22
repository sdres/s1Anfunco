"""Quantifying the peaks and troughts across layers from the center surround analysis"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'

digits = ['D2', 'D3', 'D4']

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

# Load previously extracted data
data = pd.read_csv(f'results/distanceFromPeak.csv', sep=',')

# Set styles for plotting
plt.style.use('dark_background')

paletteDigits = {'D2': '#ff180f',
                 'D3': '#00fb3b',
                 'D4': '#fffc54'}
TR = 1.9295

colors = {'#003f5c': ['#003f5c', '#8699aa', '#ffffff'],
          '#374c80': ['#374c80', '#9ba1be', '#ffffff'],
          '#7a5195': ['#7a5195', '#bca5c9', '#ffffff'],
          '#bc5090': ['#bc5090', '#e1a9c6', '#ffffff'],
          '#ef5675': ['#ef5675', '#ffafb7', '#ffffff'],
          '#ff764a': ['#ff764a', '#ffbca2', '#ffffff'],
          '#ffa600': ['#ffa600', '#ffd291', '#ffffff']
          }
hexes = ['#003f5c', '#374c80', '#7a5195', '#bc5090', '#ef5675', '#ff764a', '#ffa600']


# ===========================================================================
# Get maximum / minimum layer-specific peak and trough values (raw)

distanceList = []
layerList = []
valList = []
typeList = []
modalityList = []
timePointList = []

for modality in ['BOLD', 'VASO']:
    for layer in data['layer'].unique():
        for dist in range(2, 7):
            if modality == 'BOLD':
                tmp = data.loc[(data['modality'] == modality) &
                               (data['distVal'] == dist) &
                               (data['distFrom'] == (data['stim'])) &
                               (data['layer'] == layer) &
                               (data['volume'] >= 13) &
                               (data['volume'] <= int(13 + (30 / TR)))]
            if modality == 'VASO':
                tmp = data.loc[(data['modality'] == modality) &
                               (data['distVal'] == dist) &
                               (data['distFrom'] == (data['stim'])) &
                               (data['layer'] == layer) &
                               (data['volume'] >= 13) &
                               (data['volume'] <= int(13 + (30/TR))) &
                               (data['subject'] != 'sub-05') &
                               (data['subject'] != 'sub-06') &
                               (data['subject'] != 'sub-07') &
                               (data['subject'] != 'sub-09') &
                               (data['subject'] != 'sub-10') &
                               (data['subject'] != 'sub-18')
                ]

            maxTp = 0
            maxVal = -100

            for tp in tmp['volume'].unique():
                tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                if tmpVal > maxVal:
                    maxVal = tmpVal
                    maxTp = tp

            minTp = 0
            minVal = 100

            for tp in tmp['volume'].unique():
                tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                if tp <= maxTp:
                    continue

                if tmpVal < minVal:
                    minVal = tmpVal
                    minTp = tp

            distanceList.append(dist)
            typeList.append("Peak")
            valList.append(maxVal)
            layerList.append(layer)
            modalityList.append(modality)
            timePointList.append(maxTp)

            distanceList.append(dist)
            typeList.append("Trough")
            valList.append(minVal)
            layerList.append(layer)
            modalityList.append(modality)
            timePointList.append(minTp)


peaksData = pd.DataFrame({'data': valList,
                          'distance': distanceList,
                          'layer': layerList,
                          'type': typeList,
                          'modality': modalityList,
                          'timepoint': timePointList
                          })


# Plot peak/ trough timepoints
for k, color in enumerate(hexes):

    for j, modality in enumerate(['BOLD', 'VASO']):
        fig, axs = plt.subplots(1, 2, figsize=(9, 3))

        for i, dataType in enumerate(['Peak', 'Trough']):
            tmp = peaksData.loc[(peaksData['type'] == dataType) & (peaksData['modality'] == modality)]

            tmp['timepoint'] -= 13
            tmp['timepoint'] *= TR

            sns.barplot(tmp, ax=axs[i], x='distance', y='timepoint', hue='layer', palette=colors[color])

            axs[i].get_legend().remove()

            axs[i].set_title(dataType, fontsize=14)

            if i == 0:
                axs[i].set_ylabel(f'{modality}\nTime to Peak/Trough [s]', fontsize=14)
                if modality == 'BOLD':
                    axs[i].set_yticks(np.arange(0, 15, 2), np.arange(0, 15, 2), fontsize=12)
                if modality == 'VASO':
                    axs[i].set_yticks(np.arange(0, 28, 3), np.arange(0, 28, 3), fontsize=12)
            if i != 0:
                axs[i].set_ylabel(f'')
                if modality == 'BOLD':
                    axs[i].set_yticks(np.arange(0, 15, 2), ['']*len(np.arange(0, 15, 2)), fontsize=12)
                if modality == 'VASO':
                    axs[i].set_yticks(np.arange(0, 28, 3), ['']*len(np.arange(0, 28, 3)), fontsize=12)

            axs[i].set_xticks(range(0, 5), ['4-6', '6-8', '8-10', '10-12', '12-14'], fontsize=12)
            axs[i].set_xlabel('Distance from peak voxel [mm]', fontsize=14)

        plt.tight_layout()
        plt.savefig(f'results/group_{modality}_peakTrough_times_{k}.png', bbox_inches="tight")
        plt.show()

# Plot peak/ trough signal changes
for k, color in enumerate(hexes):
    for j, modality in enumerate(['BOLD', 'VASO']):
        fig, axs = plt.subplots(1, 2, figsize=(9, 3))
        for i, dataType in enumerate(['Peak', 'Trough']):
            tmp = peaksData.loc[(peaksData['type'] == dataType) & (peaksData['modality'] == modality)]
            sns.barplot(tmp, ax=axs[i], x='distance', y='data', hue='layer', palette=colors[color])

            axs[i].get_legend().remove()

            axs[i].set_title(dataType, fontsize=14)

            if i == 0:
                axs[i].set_ylabel(f'{modality}\nSignal change [%]', fontsize=14)
                axs[i].set_yticks(np.arange(0, 0.7, 0.1).round(decimals=1), np.arange(0, 0.7, 0.1).round(decimals=1),
                                  fontsize=12)

            if i != 0:
                axs[i].set_ylabel(f'')
                axs[i].set_yticks(np.arange(-0.6, 0.1, 0.1).round(decimals=1),
                                  np.arange(-0.6, 0.1, 0.1).round(decimals=1), fontsize=12)

            axs[i].set_xticks(range(0, 5), ['', '', '', '', ''], fontsize=12)
            axs[i].set_xlabel('', fontsize=14)
            axs[i].set_xticks(range(0, 5), ['4-6', '6-8', '8-10', '10-12', '12-14'], fontsize=12)
            axs[i].set_xlabel('Distance from peak voxel [mm]', fontsize=14)

            plt.tight_layout()
            plt.savefig(f'results/group_{modality}_peakTrough_absolute_{k}.png', bbox_inches="tight")
            plt.show()

# ===========================================================================
# For individual subs

distanceList = []
layerList = []
valList = []
typeList = []
modalityList = []
timePointList = []
subList = []


for sub in subs:
    for modality in ['BOLD', 'VASO']:
        for layer in data['layer'].unique():
            for dist in range(2, 7):
                if sub == 'sub-05':
                    continue
                if sub == 'sub-06':
                    continue
                if sub == 'sub-07':
                    continue
                if sub == 'sub-09':
                    continue
                if sub == 'sub-10':
                    continue
                if sub == 'sub-18':
                    continue

                tmp = data.loc[(data['subject'] == sub) &
                               (data['modality'] == modality) &
                               (data['distVal'] == dist) &
                               (data['distFrom'] == (data['stim'])) &
                               (data['layer'] == layer) &
                               (data['volume'] >= 13) &
                               (data['volume'] <= int(13 + (30 / TR)))]

                maxTp = 0
                maxVal = -100

                for tp in tmp['volume'].unique():
                    tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                    if tmpVal > maxVal:
                        maxVal = tmpVal
                        maxTp = tp

                minTp = 0
                minVal = 100

                for tp in tmp['volume'].unique():
                    tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                    if tp <= maxTp:
                        continue

                    if tmpVal < minVal:
                        minVal = tmpVal
                        minTp = tp

                distanceList.append(dist)
                typeList.append("Peak")
                valList.append(maxVal)
                layerList.append(layer)
                modalityList.append(modality)
                timePointList.append(maxTp)
                subList.append(sub)

                distanceList.append(dist)
                typeList.append("Trough")
                valList.append(minVal)
                layerList.append(layer)
                modalityList.append(modality)
                timePointList.append(minTp)
                subList.append(sub)


peaksData = pd.DataFrame({'subject': subList,
                         'data': valList,
                          'distance': distanceList,
                          'layer': layerList,
                          'type': typeList,
                          'modality': modalityList,
                          'timepoint': timePointList
                          })


peaksData.loc[(peaksData['type']=='Trough') & (peaksData['data']>0)]

# Plot peak/ trough signal changes
for k, color in enumerate(hexes[:1]):
    for j, modality in enumerate(['BOLD', 'VASO']):
        fig, axs = plt.subplots(1, 2, figsize=(9, 3))
        for i, dataType in enumerate(['Peak', 'Trough']):
            tmp = peaksData.loc[(peaksData['type'] == dataType) & (peaksData['modality'] == modality)]
            sns.barplot(tmp, ax=axs[i], x='distance', y='data', hue='layer', palette=colors[color])

            axs[i].get_legend().remove()

            axs[i].set_title(dataType, fontsize=14)

            if i == 0:
                axs[i].set_ylabel(f'{modality}\nSignal change [%]', fontsize=14)
                axs[i].set_yticks(np.arange(0, 0.7, 0.1).round(decimals=1), np.arange(0, 0.7, 0.1).round(decimals=1),
                                  fontsize=12)

            if i != 0:
                axs[i].set_ylabel(f'')
                axs[i].set_yticks(np.arange(-0.6, 0.1, 0.1).round(decimals=1),
                                  np.arange(-0.6, 0.1, 0.1).round(decimals=1), fontsize=12)

            axs[i].set_xticks(range(0, 5), ['', '', '', '', ''], fontsize=12)
            axs[i].set_xlabel('', fontsize=14)
            axs[i].set_xticks(range(0, 5), ['4-6', '6-8', '8-10', '10-12', '12-14'], fontsize=12)
            axs[i].set_xlabel('Distance from peak voxel [mm]', fontsize=14)

            plt.tight_layout()
            plt.savefig(f'results/group_withsubs_{modality}_peakTrough_absolute_{k}.png', bbox_inches="tight")
            plt.show()


# Investigate post stim peak

distanceList = []
layerList = []
valList = []
modalityList = []
timePointList = []

for modality in ['BOLD', 'VASO']:
    for layer in data['layer'].unique():
        for dist in range(3, 7):
            if modality == 'BOLD':
                tmp = data.loc[(data['modality'] == modality) &
                               (data['distVal'] == dist) &
                               (data['distFrom'] == (data['stim'])) &
                               (data['layer'] == layer) &
                               (data['volume'] <= 36) &
                               (data['volume'] > int(13 + (30 / TR)))]
            if modality == 'VASO':
                tmp = data.loc[(data['modality'] == modality) &
                               (data['distVal'] == dist) &
                               (data['distFrom'] == (data['stim'])) &
                               (data['layer'] == layer) &
                               (data['volume'] > int(13 + (30 / TR))) &
                               (data['volume'] <= 36) &

                               (data['subject'] != 'sub-05') &
                               (data['subject'] != 'sub-06') &
                               (data['subject'] != 'sub-07') &
                               (data['subject'] != 'sub-09') &
                               (data['subject'] != 'sub-10') &
                               (data['subject'] != 'sub-18')
                               ]

            maxTp = 0
            maxVal = -100

            for tp in tmp['volume'].unique():
                tmpVal = np.mean(tmp.loc[tmp['volume'] == tp]['data'])

                if tmpVal > maxVal:
                    maxVal = tmpVal
                    maxTp = tp

            distanceList.append(dist)
            valList.append(maxVal)
            layerList.append(layer)
            modalityList.append(modality)
            timePointList.append(maxTp)


postPeakData = pd.DataFrame({'data': valList,
                          'distance': distanceList,
                          'layer': layerList,
                          'modality': modalityList,
                          'timepoint': timePointList
                          })



# Plot peaktimepoints
for k, color in enumerate(hexes):

    for j, modality in enumerate(['BOLD', 'VASO']):
        fig, axs = plt.subplots(1, 2, figsize=(9, 3))

        tmp = postPeakData.loc[postPeakData['modality'] == modality]

        tmp['timepoint'] -= int(13 + (30 / TR))
        tmp['timepoint'] *= TR

        for i, dataType in enumerate(['timepoint', 'data']):

            sns.barplot(tmp, ax=axs[i], x='distance', y=dataType, hue='layer', palette=colors[color])
            axs[i].get_legend().remove()

            if dataType == 'timepoint':
                axs[i].set_ylabel(f'{modality}\nTime to Peak[s]', fontsize=14)
                axs[i].set_yticks(np.arange(0, 15, 2), np.arange(0, 15, 2), fontsize=12)

            if dataType == 'data':
                axs[i].set_ylabel(f'Signal change [%]', fontsize=14)
                if modality == 'BOLD':
                    axs[i].set_yticks(np.arange(0, 0.4, 0.05).round(decimals=1),
                                      np.arange(0, 0.4, 0.05).round(decimals=1), fontsize=12)
                if modality == 'VASO':
                    axs[i].set_yticks(np.arange(0, 0.3, 0.05).round(decimals=1),
                                      np.arange(0, 0.3, 0.05).round(decimals=1), fontsize=12)

            axs[i].set_xticks(range(0, 4), ['6-8', '8-10', '10-12', '12-14'], fontsize=12)
            axs[i].set_xlabel('Distance from peak voxel [mm]', fontsize=14)
        plt.suptitle('Post stimulus peak', fontsize=18)
        plt.tight_layout()
        plt.savefig(f'results/group_{modality}_poststimpeak_{k}.png', bbox_inches="tight")
        plt.show()