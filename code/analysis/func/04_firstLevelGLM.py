"""

Prepare and run first level GLM in FSL FEAT.
Assumes template files for each type of run.

"""

import os
import subprocess
import glob
import nibabel as nb
import time
import numpy as np

# Set some folder names
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

subs = ['sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for sub in subs:
    # Find all runs of participant
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-*/func/{sub}_ses-0*_task-stim_run-00*_cbv.nii.gz'))

    # Find sessions
    sessions = []
    for run in runs:
        for i in range(1, 4):
            tmp = f'ses-{i:02d}'
            if tmp in run:
                sessions.append(tmp)

    sessions = set(sessions)

    # Loop over sessions
    for ses in sessions:
        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-stim_run-00*_cbv.nii.gz'))

        for run in runs:
            base = os.path.basename(run).rsplit('.', 2)[0][:-4]

            for modality in ['BOLD', 'VASO']:

                actualData = f'{ROOT}/derivatives/{sub}/func/{base}_{modality}.nii.gz'
                runData = nb.load(actualData)

                nrVolumes = str(runData.header['dim'][4])
                print(nrVolumes)

                nrVoxels = str(np.prod(runData.header['dim'][1:5]))
                print(nrVoxels)

                replacements = {'SUBID': f'{sub}',
                                'BASE': f'{base}',
                                'ROOT': ROOT,
                                'NUMVOX': nrVoxels,
                                'NUMVOL': nrVolumes
                                }

                with open(f"{ROOT}/derivatives/designFiles/designTemplate_stim_{modality}.fsf") as infile:
                    with open(f"{ROOT}/derivatives/designFiles/fsfs/{base}_{modality}.fsf", 'w') as outfile:
                        for line in infile:
                            for src, target in replacements.items():
                                line = line.replace(src, target)
                            outfile.write(line)

# =====================================================================================================================
# Run GLMs
# =====================================================================================================================


subs = ['sub-15', 'sub-16', 'sub-17', 'sub-18']
subs = ['sub-14']


executed = 0  # Counter for parallel processes
for sub in subs:
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim_run-00*_cbv.nii.gz'))
    nrRuns = len(runs)
    print(f'Found {nrRuns} runs')

    for run in runs:
        # Set basename of run
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        for modality in ['BOLD', 'VASO']:  # Loop over modalities
        # for modality in ['VASO']:  # Loop over modalities

            # Check if the GLM for this run already ran
            file = f'{ROOT}/derivatives/{sub}/func/{base}_{modality}.feat/report_log.html'

            if os.path.exists(file):  # If yes, skip
                print(f'GLM for {base}_{modality} already ran')

            if not os.path.exists(file):  # If no, run
                print(f'Processing run {base}_{modality}')
                subprocess.run(f'feat {ROOT}/derivatives/designFiles/fsfs/{base}_{modality}.fsf &', shell=True)
                executed += 1  # Count parallel processes

            # Wait 20 minutes before starting to process next set of runs if 2 runs are being processed
            if executed >= 2:
                time.sleep(60*20)
                executed = 0  # Reset counter
