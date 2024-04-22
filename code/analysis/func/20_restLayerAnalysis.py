"""

Resting state analysis
Will replace scripts 18 and 19

"""

import nibabel as nb
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import math
import subprocess
import glob
import os
from nilearn.connectome import ConnectivityMeasure
import itertools


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

digits = ['D2', 'D3', 'D4']


# ==================================================================================================================
# Compute resting state connectivity across digits

combinations = list(itertools.combinations(digits, 2))


for sub in subs:
    print(f'Processing {sub}')
    for modality in ['VASO', 'BOLD']:

        # Load data
        data = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_layerTimecourses_clean.csv', delimiter=',')

        # deleted = np.delete(data, [0, 10, 11, 22, 23, 32], 1)
        # deleted.shape

        # Compute correlation matrix
        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([data])[0]

        # ============================================================
        # Save individual digit connectivity matrices

        for combination in combinations:
            print(combination)

            # Initialize empty matrix
            avgMatrix = np.zeros((11, 11))

            # Define number of layers
            layers = 11

            # Get number of sub-matrices
            nrSubMatrices = int(correlation_matrix.shape[0] / layers)

            # Set starting coordinates to loop over large matrix
            yOffSet = 0
            xOffSet = 0

            # Start counting added matrixes for averaging process
            matrixCount = 0

            # Loop over rows
            for x in range(nrSubMatrices):

                # Loop over columns within current row
                for y in range(nrSubMatrices):
                    # Select submatrix
                    subMatrix = correlation_matrix[yOffSet:yOffSet+layers, xOffSet:xOffSet+layers]

                    current_combination = [digits[x], digits[y]]

                    # Ignore diagonal
                    if not yOffSet == xOffSet:
                        # Add submatrix to average matrix
                        avgMatrix += subMatrix
                        # Add counter for added matrix
                        matrixCount += 1
                        # print(f'{digits[x]} - {digits[y]}')
                        print(f'{current_combination}')
                        print(f'{tuple(current_combination.sort()) == combination}')

                    # Print if matrix is skipped to verify that diagonal is not included
                    if yOffSet == xOffSet:
                        print(f'{digits[x]} - {digits[y]}')
                        print('skipping')

                    # Make sure to end loop if last row index is reached
                    if yOffSet+layers >= correlation_matrix.shape[0]:
                        break
                    # Make sure to end loop if last column index is reached
                    if xOffSet >= correlation_matrix.shape[0]:
                        break

                    # Move y-offset to next sub-matrix
                    yOffSet = yOffSet + layers

                # Jumping to new column after last row
                xOffSet = xOffSet + layers
                # Reset row index
                yOffSet = 0

            avgMatrix /= matrixCount

            # Remove outer layers
            avgMatrix = avgMatrix[1:-1, 1:-1]

            # Save null matrix
            np.savetxt(f'results/{sub}_{modality}_connectivity_condensed.csv', avgMatrix, delimiter=',')


# ==================================================================================================================
# Compute resting state connectivity across digits


for sub in subs:
    print(f'Processing {sub}')
    for modality in ['VASO', 'BOLD']:

        # Load data
        data = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_layerTimecourses_clean.csv', delimiter=',')

        # deleted = np.delete(data, [0, 10, 11, 22, 23, 32], 1)
        # deleted.shape

        # Compute correlation matrix
        correlation_measure = ConnectivityMeasure(kind='correlation')
        correlation_matrix = correlation_measure.fit_transform([data])[0]

        # ============================================================
        # Save individual digit connectivity and Average matrix across digits

        # Initialize empty matrix
        avgMatrix = np.zeros((11, 11))

        # Define number of layers
        layers = 11

        # Get number of sub-matrices
        nrSubMatrices = int(correlation_matrix.shape[0] / layers)

        # Set starting coordinates to loop over large matrix
        yOffSet = 0
        xOffSet = 0

        # Start counting added matrixes for averaging process
        matrixCount = 0

        # Loop over rows
        for x in range(nrSubMatrices):
            print('')
            print(f'new column (nr {x})')
            print(f'from {xOffSet} to {xOffSet + layers}')
            print('')

            # Loop over columns within current row
            for y in range(nrSubMatrices):

                print(f'from {yOffSet} to {yOffSet+layers}')
                # Select submatrix
                subMatrix = correlation_matrix[yOffSet:yOffSet+layers, xOffSet:xOffSet+layers]

                # Ignore diagonal
                if not yOffSet == xOffSet:
                    # Add submatrix to average matrix
                    avgMatrix += subMatrix
                    # Add counter for added matrix
                    matrixCount += 1
                    print(f'{digits[x]} - {digits[y]}')

                # Print if matrix is skipped to verify that diagonal is not included
                if yOffSet == xOffSet:
                    print(f'{digits[x]} - {digits[y]}')
                    print('skipping')

                # Make sure to end loop if last row index is reached
                if yOffSet+layers >= correlation_matrix.shape[0]:
                    break
                # Make sure to end loop if last column index is reached
                if xOffSet >= correlation_matrix.shape[0]:
                    break

                # Move y-offset to next sub-matrix
                yOffSet = yOffSet + layers

            # Jumping to new column after last row
            xOffSet = xOffSet + layers
            # Reset row index
            yOffSet = 0

        avgMatrix /= matrixCount

        # Remove outer layers
        avgMatrix = avgMatrix[1:-1, 1:-1]

        # Save null matrix
        np.savetxt(f'results/{sub}_{modality}_connectivity_condensed_normIndividuals.csv', avgMatrix, delimiter=',')


# ==================================================================================================================
# Compute resting state connectivity across digit pairs


for sub in subs:
    print(f'Processing {sub}')
    for modality in ['VASO', 'BOLD']:

        # Load data
        data = np.loadtxt(f'results/{sub}_{modality}_roi-BOLD_layerTimecourses_clean.csv', delimiter=',')


# ==================================================================================================================
# Plotting

plt.style.use('dark_background')

colors = {'BOLD': 'Reds_r',
          'VASO': 'Blues_r'}

locs = [0.5, 4.5, 8.5]
labels = ['Deep', 'Middle', 'Superficial']

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for i, sub in enumerate(subs):
    for modality in ['VASO', 'BOLD']:

        matrix = np.loadtxt(f'results/{sub}_{modality}_connectivity_condensed.csv', delimiter=',')

        fig, ax = plt.subplots()
        sns.heatmap(matrix, ax=ax, annot=False, cbar=True, square=True, cmap=colors[modality])
        ax.set_yticks(locs, labels, fontsize=12)
        ax.set_xticks(locs, labels, fontsize=12)
        ax.set_title(f'{sub} {modality} connectivity-matrix', fontsize=16)

        plt.tight_layout()
        plt.savefig(f'results/{sub}_{modality}_connectivity_condensed.png', bbox_inches="tight")
        plt.close()

