"""Script to extract laminar resting state activity from digit ROIs"""

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
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
# subs = ['sub-07', 'sub-14']
digits = ['D2', 'D3', 'D4']

pairs = {'D2D3': [11, 22, 0, 11],
         'D2D4': [22, 34, 0, 11],
         'D3D4': [22, 34, 11, 22]}

colors = {'BOLD': 'Reds',
          'VASO': 'Blues'}
#
# # for roiModality in ['VASO', 'BOLD']:
# for roiModality in ['BOLD']:
#
#     for modality in ['VASO', 'BOLD']:
#         overallMatrix = np.zeros((33, 33))
#
#         for sub in subs:
#             print(f'Processing {sub}')
#
#             # Load data
#             data = np.loadtxt(f'results/{sub}_{modality}_roi-{roiModality}_layerTimecourses_clean.csv', delimiter=',')
#
#             correlation_measure = ConnectivityMeasure(kind='correlation')
#             correlation_matrix = correlation_measure.fit_transform([data])[0]
#
#             # Mask out the major diagonal
#             np.fill_diagonal(correlation_matrix, 0)
#             overallMatrix += correlation_matrix
#
#             subCondensed = np.zeros((11, 11))
#
#             for p in pairs:
#                 subMatrix = correlation_matrix[pairs[p][0]:pairs[p][1], pairs[p][2]:pairs[p][3]]
#                 subCondensed += subMatrix
#                 fig, ax = plt.subplots()
#                 sns.heatmap(subMatrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#                 ax.set_title(f'{p}', fontsize=16)
#                 ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#                 ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
#                 ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#                 ax.set_ylabel(f'{p[:2]} layer', fontsize=14)
#                 plt.tight_layout()
#                 plt.savefig(f'results/{sub}_{p}_{modality}_roi-{roiModality}.png', bbox_inches="tight")
#                 # plt.show()
#                 plt.close()
#             subCondensed /= len(pairs)
#
#             fig, ax = plt.subplots()
#             sns.heatmap(subCondensed, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#             ax.set_title(f'{modality}\nmatrix collapsed across digits and participants', fontsize=16)
#             ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#             ax.set_xlabel(f'', fontsize=14)
#             ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#             ax.set_ylabel(f'', fontsize=14)
#             plt.tight_layout()
#             plt.savefig(f'results/{sub}_collapsed_{modality}_roi-{roiModality}.png')
#             # plt.show()
#             plt.close()
#
#         # Plot average matrix across participants
#         overallMatrix /= len(subs)
#         np.savetxt(f'results/group_{modality}_roi-{roiModality}_layerConnectivityMatrix.csv',
#                    overallMatrix,
#                    delimiter=','
#                    )
#
#         fig, ax = plt.subplots()
#         sns.heatmap(overallMatrix, ax=ax, annot=False, vmin=-0, vmax=0.5, cbar=True, square=True, cmap=colors[modality])
#         ax.set_title(f'{modality} connectivity across participants', fontsize=16)
#         # ax.set_xticks(np.linspace(0.5, 33.5, 33), range(1, 12), fontsize=12)
#         # ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
#         # ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         # ax.set_ylabel(f'{p[:2]} layer', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f'results/group_{modality}_roi-{roiModality}.png')
#         # plt.show()
#         plt.close()
#
#         subMatrix1 = overallMatrix[11:22, :11]
#         subMatrix2 = overallMatrix[22:34, :11]
#         subMatrix3 = overallMatrix[22:34, 11:22]
#
#         # Plot individual digit-pair matrices across participants
#         fig, ax = plt.subplots()
#         sns.heatmap(subMatrix1, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#         ax.set_title(f'{modality}\nmatrix collapsed across digits and participants', fontsize=16)
#         ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_xlabel(f'', fontsize=14)
#         ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_ylabel(f'', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f'results/group_D2D3collapsed_{modality}_roi-{roiModality}.png')
#         # plt.show()
#         plt.close()
#
#         fig, ax = plt.subplots()
#         sns.heatmap(subMatrix2, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#         ax.set_title(f'{modality}\nmatrix collapsed across digits and participants', fontsize=16)
#         ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_xlabel(f'', fontsize=14)
#         ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_ylabel(f'', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f'results/group_D2D4collapsed_{modality}_roi-{roiModality}.png')
#         # plt.show()
#         plt.close()
#
#         fig, ax = plt.subplots()
#         sns.heatmap(subMatrix3, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#         ax.set_title(f'{modality}\nmatrix collapsed across digits and participants', fontsize=16)
#         ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_xlabel(f'', fontsize=14)
#         ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_ylabel(f'', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f'results/group_D3D4collapsed_{modality}_roi-{roiModality}.png')
#         # plt.show()
#         plt.close()
#
#         # Average all digit-pair matrices across participants
#         subMatrix1 += subMatrix2
#         subMatrix1 += subMatrix3
#         subMatrix1 /= 3
#
#         fig, ax = plt.subplots()
#         sns.heatmap(subMatrix1, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
#         ax.set_title(f'{modality}\nmatrix collapsed across digits and participants', fontsize=16)
#         ax.set_xticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_xlabel(f'', fontsize=14)
#         ax.set_yticks(np.linspace(0.5, 10.5, 11), range(1, 12), fontsize=12)
#         ax.set_ylabel(f'', fontsize=14)
#         plt.tight_layout()
#         plt.savefig(f'results/group_collapsed_{modality}_roi-{roiModality}.png')
#         # plt.show()
#         plt.close()

# =================================================================================
# Including variability

roiModality = 'BOLD'
subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-14', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-07']

for modality in ['VASO', 'BOLD']:
    overallMatrix = np.zeros((33, 33, len(subs)))

    for i, sub in enumerate(subs):
        print(f'Processing {sub}')

        # Load data
        data = np.loadtxt(f'results/{sub}_{modality}_roi-{roiModality}_layerTimecourses_clean.csv', delimiter=',')
        nullData = np.loadtxt(f'results/{sub}_{modality}_nullMatrix.csv', delimiter=',')

        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([data])[0]

        # Mask out the major diagonal
        np.fill_diagonal(correlation_matrix, 0)
        overallMatrix[:, :, i] = correlation_matrix

        subCondensed = np.zeros((11, 11))
        subCondensedNorm = np.zeros((11, 11))
        subCondensedCorrected = np.zeros((11, 11))

        for p in pairs:
            subMatrix = correlation_matrix[pairs[p][0]:pairs[p][1], pairs[p][2]:pairs[p][3]]

            subCondensed += subMatrix

            minVal = np.min(subMatrix)
            maxVal = np.max(subMatrix)
            normMatrix = ((subMatrix.copy() - minVal) / (maxVal - minVal))

            subCondensedNorm += normMatrix

            fig, ax = plt.subplots()
            sns.heatmap(subMatrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
            ax.set_title(f'{sub} {p}', fontsize=16)
            ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
            ax.set_ylabel(f'{p[:2]} layer', fontsize=14)

            ax.set_yticks(locs, labels, fontsize=12)
            ax.set_xticks(locs, labels, fontsize=12)

            plt.tight_layout()
            plt.savefig(f'results/{sub}_{p}_{modality}_roi-{roiModality}.png', bbox_inches="tight")
            # plt.show()
            plt.close()

            subMatrix /= nullData
            subCondensedCorrected += subMatrix

            fig, ax = plt.subplots()
            sns.heatmap(subMatrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
            ax.set_title(f'{sub} {p}', fontsize=16)
            ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
            ax.set_ylabel(f'{p[:2]} layer', fontsize=14)

            ax.set_yticks(locs, labels, fontsize=12)
            ax.set_xticks(locs, labels, fontsize=12)

            plt.tight_layout()
            plt.savefig(f'results/{sub}_{p}_{modality}_roi-{roiModality}_correctedByDivision.png', bbox_inches="tight")
            # plt.show()
            plt.close()


            fig, ax = plt.subplots()
            sns.heatmap(normMatrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
            ax.set_title(f'{sub} {p}', fontsize=16)
            ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
            ax.set_ylabel(f'{p[:2]} layer', fontsize=14)

            ax.set_yticks(locs, labels, fontsize=12)
            ax.set_xticks(locs, labels, fontsize=12)

            plt.tight_layout()
            plt.savefig(f'results/{sub}_{p}_{modality}_roi-{roiModality}_normalised.png', bbox_inches="tight")
            # plt.show()
            plt.close()


        subCondensed /= len(pairs)
        subCondensedCorrected /= len(pairs)
        subCondensedNorm /= len(pairs)

        fig, ax = plt.subplots()
        sns.heatmap(subCondensed, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
        ax.set_title(f'{sub} collapsed across digits', fontsize=16)
        ax.set_xlabel(f'', fontsize=14)
        ax.set_ylabel(f'', fontsize=14)

        ax.set_yticks(locs, labels, fontsize=12)
        ax.set_xticks(locs, labels, fontsize=12)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_collapsed_{modality}_roi-{roiModality}.png')
        # plt.show()
        plt.close()

        fig, ax = plt.subplots()
        sns.heatmap(subCondensedCorrected, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
        ax.set_title(f'{sub} collapsed across digits (corrected)', fontsize=16)
        ax.set_xlabel(f'', fontsize=14)
        ax.set_ylabel(f'', fontsize=14)

        ax.set_yticks(locs, labels, fontsize=12)
        ax.set_xticks(locs, labels, fontsize=12)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_collapsedCorrected_{modality}_roi-{roiModality}.png')
        # plt.show()
        plt.close()

        fig, ax = plt.subplots()
        sns.heatmap(subCondensedNorm, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
        ax.set_title(f'{sub} collapsed across digits (corrected)', fontsize=16)
        ax.set_xlabel(f'', fontsize=14)
        ax.set_ylabel(f'', fontsize=14)

        ax.set_yticks(locs, labels, fontsize=12)
        ax.set_xticks(locs, labels, fontsize=12)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_collapsedNormalised_{modality}_roi-{roiModality}.png')
        # plt.show()
        plt.close()


    # Plot average matrix across participants
    averageOverall = np.mean(overallMatrix, axis=-1)

    np.savetxt(f'results/group_{modality}_roi-{roiModality}_layerConnectivityMatrix.csv',
               overallMatrix,
               delimiter=','
               )

# =================================================================================
# Plotting

roiModality = 'BOLD'

locs = [0.5, 5.5, 10.5]
labels = ['Deep', 'Middle', 'Superficial']

fig, ax = plt.subplots()
sns.heatmap(overallMatrix, ax=ax, annot=False, vmin=-0, vmax=0.5, cbar=True, square=True, cmap=colors[modality])
ax.set_title(f'connectivity across participants', fontsize=16)
# ax.set_xlabel(f'{p[2:]} layer', fontsize=14)
# ax.set_ylabel(f'{p[:2]} layer', fontsize=14)

ax.set_yticks(locs, labels, fontsize=12)
ax.set_xticks(locs, labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_{modality}_roi-{roiModality}.png')
# plt.show()
plt.close()

subMatrix1 = overallMatrix[11:22, :11]
subMatrix2 = overallMatrix[22:34, :11]
subMatrix3 = overallMatrix[22:34, 11:22]

# Plot individual digit-pair matrices across participants
fig, ax = plt.subplots()
sns.heatmap(subMatrix1, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
ax.set_title(f'D2-D3 matrix across participants', fontsize=16)
ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks(locs, labels, fontsize=12)
ax.set_xticks(locs, labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_D2D3collapsed_{modality}_roi-{roiModality}.png')
# plt.show()
plt.close()

fig, ax = plt.subplots()
sns.heatmap(subMatrix2, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
ax.set_title(f'D2-D4 matrix across participants', fontsize=16)
ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks(locs, labels, fontsize=12)
ax.set_xticks(locs, labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_D2D4collapsed_{modality}_roi-{roiModality}.png')
# plt.show()
plt.close()

fig, ax = plt.subplots()
sns.heatmap(subMatrix3, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
ax.set_title(f'D3-D4 matrix across participants', fontsize=16)
ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks(locs, labels, fontsize=12)
ax.set_xticks(locs, labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_D3D4collapsed_{modality}_roi-{roiModality}.png')
# plt.show()
plt.close()

# Average all digit-pair matrices across participants
avg = np.zeros((11, 11))
avg += subMatrix1
avg += subMatrix2
avg += subMatrix3
avg /= 3

fig, ax = plt.subplots()
sns.heatmap(avg, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
ax.set_title(f'matrix collapsed across digits and participants', fontsize=16)
ax.set_xlabel(f'', fontsize=14)
ax.set_ylabel(f'', fontsize=14)

ax.set_yticks(locs, labels, fontsize=12)
ax.set_xticks(locs, labels, fontsize=12)

plt.tight_layout()
plt.savefig(f'results/group_collapsed_{modality}_roi-{roiModality}.png')
# plt.show()
plt.close()

# Get diagonals
diag = subMatrix2.diagonal()[1:-1]
x = np.arange(len(diag))
plt.plot(x, diag, label=modality)
plt.legend()
plt.savefig(f'results/group_D2D4collapsed_{modality}_roi-{roiModality}_profile.png')
plt.show()
