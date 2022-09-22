'''

Populates .tsv files in subject-specific BIDS directory

'''


import pandas as pd
import glob
import os

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

subs = ['sub-07']

for sub in subs:
    print(sub)

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    runs = sorted(glob.glob(f'{ROOT}/{sub}/*/func/{sub}_*_task-stim*run-00*_cbv.nii.gz'))

    sessions = []

    for run in runs:
        for sesNr in range(1,4):
            currSes = f'ses-0{sesNr}'
            if currSes in run:
                sessions.append(currSes)

    sessions = set(sessions)  # Remove duplicates

    for ses in sessions:

        for run in runs:
            base = os.path.basename(run).rsplit('.', 2)[0][:-4]

            events = pd.DataFrame(columns = ['onset', 'duration', 'trial_type'])

            for digit in ['D2', 'D3', 'D4']:
                try:
                    file = f'{ROOT}/derivatives/designFiles/stimulationPatterns/{sub}/{base}_{digit}.txt'
                    tmp = pd.read_csv(file, sep=' ', names = ['onset', 'duration', 'trial_type'])

                except:
                    print('old naming convention')
                    fileNew = file[:-26] + '0' + file[-26:]
                    tmp = pd.read_csv(fileNew, sep=' ', names = ['onset', 'duration', 'trial_type'])


                tmp['trial_type'] = [digit] * 4  # 4 repetitions per finger
                events = pd.concat((events, tmp))
                events.sort_values('onset', inplace = True)

            events.to_csv(f'{ROOT}/{sub}/{ses}/func/{base}_events.tsv', sep = ' ',index=False)
