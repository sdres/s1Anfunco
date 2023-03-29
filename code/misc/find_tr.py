'''

Takes logfile as input ant reads TR from trigger times.
Assumes keypress: 5 as trigger

'''

import re
import numpy as np


logfile = '/Volumes/extData/S1ANFUNCO/Nifti/derivatives/logFiles/sub-06_ses-001_task-stim_run-001_cbv.log'

f = open(logfile)

df = pd.read_csv(f, skiprows=2, sep='\t', lineterminator='\n',
                 usecols=["Code", 'Time'])
df.rename(columns={'Time': 'start_time', 'Code': 'Type'}, inplace=True)

df["start_time"] = pd.to_numeric(df["start_time"])

first_pulse_time = df.at[1, 'start_time']
last_pulse_time = int(df.at[df.last_valid_index(), 'start_time'])

df['start_time'] = df['start_time'].sub(first_pulse_time)  # normalize to first pulse-time
df['start_time'] = df['start_time'].div(10000)  # convert to seconds

triggers = df.loc[df['Type']=='44']

triggerTimes = triggers['start_time'].to_numpy()

triggersSubtracted = []
for n in range(len(triggerTimes) - 1):
    triggersSubtracted.append(float(triggerTimes[n + 1]) - float(triggerTimes[n]))

meanFirstTriggerDur = np.mean(triggersSubtracted[::2])
meanSecondTriggerDur = np.mean(triggersSubtracted[1::2])

# find mean trigger-time
meanTriggerDur = (meanFirstTriggerDur + meanSecondTriggerDur) / 2

TRpair = meanFirstTriggerDur + meanSecondTriggerDur

inversionDelay = 0.650

inversionTime = (meanFirstTriggerDur-inversionDelay)/2 + inversionDelay

