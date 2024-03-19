"""Script to correct resting state matrices with null-matrix"""

import nibabel as nb
import numpy as np
import glob
from nilearn.connectome import ConnectivityMeasure
from nilearn import plotting
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

plt.style.use('dark_background')


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-05']
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10']
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-07', 'sub-14']

# subs = ['sub-18']

locs = [0.5, 5.5, 10.5]
labels = ['Deep', 'Middle', 'Superficial']

digits = ['D2', 'D3', 'D4']

pairs = {'D2D3': [11, 22, 0, 11],
         'D2D4': [22, 34, 0, 11],
         'D3D4': [22, 34, 11, 22]}

colors = {'BOLD': 'Reds',
          'VASO': 'Blues'}


roiModality = 'BOLD'
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-14', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-07']

# ======================================================================
# Take relative difference between null and hypothesis matrix

for modality in ['VASO', 'BOLD']:

    across = np.zeros((11, 11))
    null = np.zeros((11, 11))
    relDiffParts = np.zeros((11, 11))

    for i, sub in enumerate(subs):
        print(f'Processing {sub}')

        # Load data
        data = np.loadtxt(f'results/{sub}_{modality}_roi-{roiModality}_layerTimecourses_clean.csv', delimiter=',')
        # deleted = np.delete(data, [0, 10, 11, 22, 23, 32], 1)
        # deleted.shape
        nullData = np.loadtxt(f'results/{sub}_{modality}_nullMatrix.csv', delimiter=',')
        # nullData = nullData[1:-1, 1:-1]

        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([data])[0]

        # matrix = np.triu(correlation_matrix)

        # fig, ax = plt.subplots()
        # sns.heatmap(correlation_matrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality], mask=matrix)
        # ax.set_title(f'{modality}', fontsize=16)
        # ax.xaxis.set_tick_params(labeltop=True)
        # ax.xaxis.set_tick_params(labelbottom=False)
        # ax.xaxis.tick_top()
        #
        # # ax.set_xlabel(f'', fontsize=14)
        # # ax.set_ylabel(f'', fontsize=14)
        #
        # ax.set_yticks([4.5, 13.5, 22.5], ["D2", 'D3', 'D4'], fontsize=12)
        # ax.set_xticks([4.5, 13.5, 22.5], ["D2", 'D3', 'D4'], fontsize=12)
        # #
        # # ax.set_yticks([0.5, 4.5, 13.5, 22.5, 26.5], ['Deep', "D2", 'D3', 'D4', 'Superficial'], fontsize=12)
        # # ax.set_xticks([0.5, 4.5, 13.5, 22.5, 26.5], ['Deep', "D2", 'D3', 'D4', 'Superficial'], fontsize=12)
        #
        # plt.tight_layout()
        # plt.savefig(f'results/{sub}_fullmatrix_{modality}.png')
        # plt.show()
        # # plt.close()





        # Average digit connectivity matrices
        subCondensed = np.zeros((11, 11))
        for p in pairs:
            subMatrix = correlation_matrix[pairs[p][0]:pairs[p][1], pairs[p][2]:pairs[p][3]]
            subCondensed += subMatrix.copy()

        subCondensed /= len(pairs)
        # subCondensed = subCondensed[1:-1, 1:-1]

        # =====================================================================================
        # Compute rel Diff per participant

        # Find uppler and lower thresh
        minAcrossNorm, maxAcrossNorm = np.percentile(subCondensed, (1, 99))
        minNullNorm, maxNullNorm = np.percentile(nullData, (1, 99))

        # Clip to extremes
        subCondensed[subCondensed < minAcrossNorm] = minAcrossNorm
        subCondensed[subCondensed > maxAcrossNorm] = maxAcrossNorm

        nullData[nullData < minNullNorm] = minNullNorm
        nullData[nullData > maxNullNorm] = maxNullNorm

        # Normalise to 0-1 range
        acrossNorm = (subCondensed - minAcrossNorm) / (maxAcrossNorm - minAcrossNorm)
        nullNorm = (nullData - minNullNorm) / (maxNullNorm - minNullNorm)

        relDiffNormSub = (acrossNorm - nullNorm) / (acrossNorm + nullNorm)

        fig, ax = plt.subplots()
        sns.heatmap(relDiffNormSub[1:-1, 1:-1], ax=ax, annot=False, cbar=True, vmax=0.4, vmin=-0.4, square=True, cmap='bwr')
        ax.set_title(f'{sub} {modality}', fontsize=16)

        ax.xaxis.set_tick_params(labeltop=True)
        ax.xaxis.set_tick_params(labelbottom=False)
        ax.xaxis.tick_top()

        ax.set_xlabel(f'', fontsize=14)
        ax.set_ylabel(f'', fontsize=14)

        ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
        ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_relDiffNormForSubs_{modality}.png')
        plt.show()
        # plt.close()



        relDiffParts += relDiffNormSub
        across += subCondensed
        null += nullData

    across /= len(subs)
    null /= len(subs)
    relDiffParts /= len(subs)

    # across matrix without processing
    fig, ax = plt.subplots()
    sns.heatmap(across[1:-1, 1:-1], ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
    ax.set_title(f'{modality}', fontsize=16)

    ax.xaxis.set_tick_params(labeltop=True)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.xaxis.tick_top()

    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
    ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_acrossMat_{modality}.png')
    plt.show()
    # plt.close()




    # Plot matrix for subject-wise relative differences
    fig, ax = plt.subplots()
    sns.heatmap(relDiffParts, ax=ax, annot=False, cbar=True, vmax=0.4, vmin=-0.4, square=True, cmap='bwr')
    ax.set_title(f'{modality}', fontsize=16)

    ax.xaxis.set_tick_params(labeltop=True)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.xaxis.tick_top()

    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
    ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_relDiffNormForSubs_{modality}.png')
    plt.show()
    # plt.close()

    # Find uppler and lower thresh
    minAcrossNorm, maxAcrossNorm = np.percentile(across, (1, 99))
    minNullNorm, maxNullNorm = np.percentile(null, (1, 99))

    # Clip to extremes
    across[across < minAcrossNorm] = minAcrossNorm
    across[across > maxAcrossNorm] = maxAcrossNorm

    null[null < minNullNorm] = minNullNorm
    null[null > maxNullNorm] = maxNullNorm

    # Normalise to 0-1 range
    acrossNorm = (across - minAcrossNorm) / (maxAcrossNorm - minAcrossNorm)
    nullNorm = (null - minNullNorm) / (maxNullNorm - minNullNorm)

    relDiffNorm = (acrossNorm - nullNorm) / (acrossNorm + nullNorm)
    relDiff = (across - null) / (across + null)

    fig, ax = plt.subplots()
    sns.heatmap(relDiffNorm, ax=ax, annot=False, cbar=True, vmax=0.5, vmin=-0.5, square=True, cmap='bwr')
    ax.set_title(f'{modality}', fontsize=16)
    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks(locs, labels, fontsize=12)
    ax.set_xticks(locs, labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_relDiffNorm_{modality}.png')
    plt.show()
    # plt.close()

    fig, ax = plt.subplots()
    sns.heatmap(relDiff, ax=ax, annot=False, cbar=True, vmax=1, vmin=0, square=True, cmap='bwr')
    ax.set_title(f'{modality}', fontsize=16)
    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks(locs, labels, fontsize=12)
    ax.set_xticks(locs, labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_relDiff_{modality}.png')
    plt.show()
    # plt.close()

    # Remove outer layers

    relDiffNorm_minOuter = relDiffNorm[1:-1, 1:-1]

    fig, ax = plt.subplots()
    sns.heatmap(relDiffNorm_minOuter, ax=ax, annot=False, cbar=True, vmax=0.4, vmin=-0.4, square=True, cmap='bwr')
    ax.set_title(f'{modality}', fontsize=16)

    ax.xaxis.set_tick_params(labeltop=True)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.xaxis.tick_top()

    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
    ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_relDiffNorm_{modality}_removeOuter.png')
    plt.show()
    # plt.close()

    null_minOuter = null[1:-1, 1:-1]

    fig, ax = plt.subplots()
    sns.heatmap(null_minOuter, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
    ax.set_title(f'{modality}', fontsize=16)

    ax.xaxis.set_tick_params(labeltop=True)
    ax.xaxis.set_tick_params(labelbottom=False)
    ax.xaxis.tick_top()

    ax.set_xlabel(f'', fontsize=14)
    ax.set_ylabel(f'', fontsize=14)

    ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
    ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

    plt.tight_layout()
    plt.savefig(f'results/group_null_{modality}_removeOuter.png')
    plt.show()
    # plt.close()

# ======================================================================
# Compare BOLD and VASO null matrix


null_bold = np.zeros((9, 9))
null_vaso = np.zeros((9, 9))
null_bold_norm = np.zeros((9, 9))
null_vaso_norm = np.zeros((9, 9))

for i, sub in enumerate(subs):
    print(f'Processing {sub}')

    for modality in ['VASO', 'BOLD']:

        # Load data
        null = np.loadtxt(f'results/{sub}_{modality}_nullMatrix.csv', delimiter=',')
        # Remove outer layers
        null = null[1: -1, 1: -1]

        if modality == 'BOLD':
            null_bold += null
        if modality == 'VASO':
            null_vaso += null

        # Find upper and lower thresh
        minNullNorm, maxNullNorm = np.percentile(null, (1, 99))

        # Clip to extremes
        null[null < minNullNorm] = minNullNorm
        null[null > maxNullNorm] = maxNullNorm

        # Normalise to 0-1 range
        nullNorm = (null - minNullNorm) / (maxNullNorm - minNullNorm)

        if modality == 'BOLD':
            null_bold_norm += nullNorm
        if modality == 'VASO':
            null_vaso_norm += nullNorm


null_bold /= len(subs)
null_vaso /= len(subs)

null_bold_norm /= len(subs)
null_vaso_norm /= len(subs)

# Compute relative difference between bold and vaso null matrix
relDiff = (null_bold - null_vaso) / (null_bold + null_vaso)

fig, ax = plt.subplots()
sns.heatmap(relDiff, ax=ax, annot=False, vmax=0.4, vmin=-0.4, cbar=True, square=True, cmap='bwr')
ax.set_title(f'BOLD null > VASO null', fontsize=16)

ax.xaxis.set_tick_params(labeltop=True)
ax.xaxis.set_tick_params(labelbottom=False)
ax.xaxis.tick_top()

ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_relDiffboldvasonull_{modality}.png')
plt.show()


# Compute relative difference between bold and vaso null matrix for normalised
relDiff = (null_bold_norm - null_vaso_norm) / (null_bold_norm + null_vaso_norm)

fig, ax = plt.subplots()
sns.heatmap(relDiff, ax=ax, annot=False, vmax=0.4, vmin=-0.4, cbar=True, square=True, cmap='bwr')
ax.set_title(f'BOLD null > VASO null normalised', fontsize=16)

ax.xaxis.set_tick_params(labeltop=True)
ax.xaxis.set_tick_params(labelbottom=False)
ax.xaxis.tick_top()

ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_relDiffboldvasonull_{modality}_norm.png')
plt.show()


# Normalise after

# Find upper and lower thresh
minNullNorm, maxNullNorm = np.percentile(null_bold, (1, 99))

# Clip to extremes
null_bold[null_bold < minNullNorm] = minNullNorm
null_bold[null_bold > maxNullNorm] = maxNullNorm

# Normalise to 0-1 range
nullNormBOLD = (null_bold - minNullNorm) / (maxNullNorm - minNullNorm)

# Find upper and lower thresh
minNullNorm, maxNullNorm = np.percentile(null_vaso, (1, 99))

# Clip to extremes
null_vaso[null_vaso < minNullNorm] = minNullNorm
null_vaso[null_vaso > maxNullNorm] = maxNullNorm

# Normalise to 0-1 range
nullNormVASO = (null_vaso - minNullNorm) / (maxNullNorm - minNullNorm)

relDiff = (nullNormBOLD - nullNormVASO) / (nullNormBOLD + nullNormVASO)

fig, ax = plt.subplots()
sns.heatmap(relDiff, ax=ax, annot=False, vmax=0.4, vmin=-0.4, cbar=True, square=True, cmap='bwr')
ax.set_title(f'BOLD null > VASO null normalised', fontsize=16)

ax.xaxis.set_tick_params(labeltop=True)
ax.xaxis.set_tick_params(labelbottom=False)
ax.xaxis.tick_top()

ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks([0.5, 4.5, 8.5], labels, fontsize=12)
ax.set_xticks([0.5, 4.5, 8.5], labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_relDiffboldvasonull_{modality}_normLast.png')
plt.show()