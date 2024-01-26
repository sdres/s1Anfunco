"""

Collect statistical maps of all participants in a folder

"""

import glob
import shutil


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'


subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']


contrasts = {1: 'D2vsAll',
             2: 'D3vsAll',
             3: 'D4vsAll',
             4: 'D2vsRest',
             5: 'D3vsRest',
             6: 'D4vsRest'
             }

for sub in subs:
    print(sub)
    # Look for all runs in all sessions
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-*/func/{sub}_ses-*_task-*run-00*_cbv.nii.gz'))

    # Look for sessions of this participant
    sessions = []
    for run in runs:
        for i in range(1, 4):
            if f'ses-0{i}' in run:
                ses = f'ses-0{i}'
                sessions.append(ses)

    sessions = set(sessions)  # Remove duplicates

    for ses in sessions:
        # Look for runs of current session
        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-stim_run-00*_cbv.nii.gz'))

        statsDir = f'{ROOT}/derivatives/{sub}/func/statMaps'

        if not os.path.exists(statsDir):
            os.makedirs(statsDir)
            print("Stats directory is created")

        nrRuns = len(runs)

        funcDir = f'{ROOT}/derivatives/{sub}/func'

        for modality in ['BOLD', 'VASO']:
            if nrRuns == 1:  # Get data from first level feat folder
                for contrastNr in contrasts:
                    statMap = f'{funcDir}/{sub}_{ses}_task-stim_run-001_{modality}.feat/stats/zstat{contrastNr}.nii.gz'
                    outName = f'{statsDir}/{sub}_{contrasts[contrastNr]}_{modality}.nii.gz'
                    shutil.copyfile(statMap, outName)
                    # print(statMap)
                    # print(outName)

            elif nrRuns >= 2:
                for contrastNr in contrasts:
                    statMap = f'{funcDir}/{sub}_{ses}_task-stim_run-secondLevel_{modality}.gfeat/cope{contrastNr}.feat/stats/zstat1.nii.gz'
                    outName = f'{statsDir}/{sub}_{contrasts[contrastNr]}_{modality}.nii.gz'
                    shutil.copyfile(statMap, outName)
                    # print(statMap)
                    # print(outName)


# ===========================================
# For sub-12
# ===========================================

subs = ['sub-12']


contrasts = {1: 'D1vsAll',
             2: 'D2vsAll',
             3: 'D3vsAll',
             4: 'D4vsAll',
             5: 'D5vsAll',
             6: 'D1vsRest',
             7: 'D2vsRest',
             8: 'D3vsRest',
             9: 'D4vsRest',
             10: 'D5vsRest'
             }

for sub in subs:
    print(sub)
    # Look for all runs in all sessions
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-*/func/{sub}_ses-*_task-*run-00*_cbv.nii.gz'))

    # Look for sessions of this participant
    sessions = []
    for run in runs:
        for i in range(1, 4):
            if f'ses-0{i}' in run:
                ses = f'ses-0{i}'
                sessions.append(ses)

    sessions = set(sessions)  # Remove duplicates

    for ses in sessions:
        # Look for runs of current session
        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-stim_run-00*_cbv.nii.gz'))

        statsDir = f'{ROOT}/derivatives/{sub}/func/statMaps'

        if not os.path.exists(statsDir):
            os.makedirs(statsDir)
            print("Stats directory is created")

        nrRuns = len(runs)

        funcDir = f'{ROOT}/derivatives/{sub}/func'

        for modality in ['BOLD', 'VASO']:
            if nrRuns == 1:  # Get data from first level feat folder
                for contrastNr in contrasts:
                    statMap = f'{funcDir}/{sub}_{ses}_task-stim_run-001_{modality}.feat/stats/zstat{contrastNr}.nii.gz'
                    outName = f'{statsDir}/{sub}_{contrasts[contrastNr]}_{modality}.nii.gz'
                    shutil.copyfile(statMap, outName)
                    # print(statMap)
                    # print(outName)

            elif nrRuns >= 2:
                for contrastNr in contrasts:
                    statMap = f'{funcDir}/{sub}_{ses}_task-stim_run-secondLevel_{modality}.gfeat/cope{contrastNr}.feat/stats/zstat1.nii.gz'
                    outName = f'{statsDir}/{sub}_{contrasts[contrastNr]}_{modality}.nii.gz'
                    shutil.copyfile(statMap, outName)
                    # print(statMap)
                    # print(outName)

