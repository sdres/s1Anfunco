"""
Plotting motion traces for visual inspection.

"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import numpy as np

subs = ['sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

# subs = ['sub-15']

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

v1Palette = {
    'notnulled': 'tab:orange',
    'nulled': 'tab:blue'}


plt.style.use('dark_background')


for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim_run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')
        motionTraces = []
        motionNames = []
        volumes = []
        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'

        for modality in ['nulled', 'notnulled']:

            data = pd.read_csv(os.path.join(motionDir, f'{base}_{modality}_motionRegressors.csv'))

            for col in data.columns:
                tmp = data[col].to_numpy()

                for i, val in enumerate(tmp, start=1):
                    motionTraces.append(val)
                    motionNames.append(f'{col} {modality}')
                    volumes.append(i)

        data_dict = {
            'volume': volumes,
            'Motion': motionTraces,
            'name': motionNames,
        }

        pd_ses = pd.DataFrame(data=data_dict)
        pd_ses.to_csv(os.path.join(motionDir, f'{base}_motionSummary.csv'), index=False)


width = 2  # Set linewidth

for sub in subs:

    outDir = os.path.join(os.getcwd(), f'results/motion/{sub}')

    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("Output directory is created")

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'

        data = pd.read_csv(os.path.join(motionDir, f'{base}_motionSummary.csv'))

        fig, axes = plt.subplots(1, 2, sharex=True, figsize=(30, 6))
        plt.suptitle(f'{base} Motion Summary', fontsize=24)

        for i, type in enumerate(['T', 'R']):

            if type == 'T':
                legend = False
            if type == 'R':
                legend = True

            for modality, cmap in zip(['nulled', 'notnulled'], ['Set1', 'Set2']):

                if modality == 'nulled':
                    tmpData = data.loc[(data['name'].str.contains(type)) & (-data['name'].str.contains('notnulled'))]
                else:
                    tmpData = data.loc[(data['name'].str.contains(type)) & (data['name'].str.contains('notnulled'))]

                sns.lineplot(ax=axes[i],
                             x='volume',
                             y='Motion',
                             data=tmpData,
                             hue='name',
                             palette=cmap,
                             linewidth=width,
                             legend=legend
                             )

        axes[0].set_ylabel("Translation [mm]", fontsize=24)
        axes[1].set_ylabel("Rotation [radians]", fontsize=24)

        # Colors are the same for both plots so we only need one legend
        axes[1].legend(fontsize=20, loc='center left', bbox_to_anchor=(1, 0.5))

        for j in range(2):
            axes[j].tick_params(axis='both', labelsize=20)
            axes[j].set_xlabel("Volume", fontsize=24)

        plt.savefig(f'{outDir}/{base}_motion.jpg', bbox_inches = 'tight', pad_inches = 0)
        plt.show()

# ============================================
# Plot FDs

for sub in subs:

    outDir = os.path.join(os.getcwd(), f'results/motion/{sub}')

    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("Output directory is created")

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim_run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'

        FDs = pd.read_csv(os.path.join(motionDir, f'{base}_FDs.csv'))

        fig = plt.figure(figsize=(20, 5))
        sns.lineplot(data=FDs, x='volume', y='FD', hue='modality', linewidth=width, palette=v1Palette)

        if np.max(FDs['FD']) < 1:
            plt.ylim(0, 1)

        plt.ylabel('FD [mm]', fontsize=24)
        plt.xlabel('Volume', fontsize=24)
        plt.title(f"{base}", fontsize=24, pad=20)
        plt.axhline(0.75, color='gray', linestyle='--', label='In-plane res.')

        plt.legend(fontsize=20, loc='upper right')

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.savefig(f'{outDir}/{base}_FDs.png', bbox_inches='tight', pad_inches=0)
        plt.show()


# =============================================
# Plot number of motion volumes per participant

FDall = pd.DataFrame()

runList = []
for sub in subs:

    outDir = os.path.join(os.getcwd(), f'results/motion/{sub}')

    if not os.path.exists(outDir):
        os.makedirs(outDir)
        print("Output directory is created")

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-stim*run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'

        FDs = pd.read_csv(os.path.join(motionDir, f'{base}_FDs.csv'))
        FDall = pd.concat([FDall, FDs])

        for i in range(len(FDs)):
            runList.append(base)


FDall['run'] = runList
nrVolsTotal = len(FDall)
print(f'the total numer of volumes is {nrVolsTotal}')

PALETTE = {'notnulled': 'tab:orange', 'nulled': 'tab:blue'}

nrVolsExcessive = 0
for modality in ['nulled', 'notnulled']:
    # Select volumes with excessive motion
    tmp = FDall.loc[(FDall['FD'] > 0.75) & (FDall['modality'] == modality)]

    maxMotion = np.max(tmp['FD'])
    print(f'maximum FD for {modality} is {maxMotion}')

    nrVolsExcessive += len(tmp)
    print(f'For {modality}, we found {len(tmp)} volumes with motion > 0.75mm')

    nrRunsWithExcessive = len(tmp['run'].unique())

    print(f'For {modality}, {nrRunsWithExcessive} runs were affected')

percentVols = nrVolsExcessive / nrVolsTotal * 100
print(f'{percentVols:2f}% of volumes show motion above voxel size')

tmp = FDall.loc[(FDall['FD'] > 0.75)]

fig, (ax1) = plt.subplots(1, 1, figsize=(7.5, 5))
# sns.histplot(tmp, x='FD', linewidth=1, color='tab:red', bins=20)
sns.histplot(tmp, x='FD', linewidth=1, hue='modality', multiple='dodge', bins=20, palette=PALETTE)
ax1.set_ylabel(r'# Volumes', fontsize=24)
ax1.set_xlabel(r'FD [mm]', fontsize=24)
ax1.yaxis.set_tick_params(labelsize=18)
ax1.xaxis.set_tick_params(labelsize=18)

plt.legend(title='Modality', loc='upper right', labels=['Notnulled', 'Nulled'], fontsize=18, title_fontsize=20)

fig.tight_layout()
plt.savefig(f'./results/runsGreater075mmFD.png', bbox_inches="tight")
plt.show()


# Check how many volumes were scrubbed due to motion outliers

scrubbedList = []
runList = []
modalityList = []
overallScrubbed = 0
for sub in subs:
    motionFolder = f'{ROOT}/derivatives/{sub}/func/motionParameters'
    for modality in ['nulled', 'notnulled']:
        files = sorted(glob.glob(f'{motionFolder}/*stim*/sub-*_ses-0*_task-stim_run-00*_{modality}_motionOutliers.txt'))
        for file in files:
            base = os.path.basename(file).rsplit('.', 2)[0][:-15]

            data = np.loadtxt(file)
            nrVols = np.sum(data)

            overallScrubbed += nrVols
            runList.append(base)
            scrubbedList.append(nrVols)
            modalityList.append(modality)
            print(f'{base} {nrVols}')


scrubbedData = pd.DataFrame({'run': runList, 'scrubbed': scrubbedList, 'modality': modalityList})

PALETTE = {'notnulled': 'tab:orange', 'nulled': 'tab:blue'}

fig, ax = plt.subplots(1, 1, figsize=(7.5, 5))
sns.histplot(data=scrubbedData, x='scrubbed', hue='modality', multiple='dodge', linewidth=1, palette=PALETTE)

plt.ylabel('# Runs', fontsize=24)
plt.yticks(np.arange(0, 8), fontsize=18)

ticks = np.arange(4, 20, 2)
plt.xticks(ticks, fontsize=18)
# plt.xlim([0, 50])

plt.xlabel(f'# Censored volumes', fontsize=18)

plt.legend(title='Modality', loc='upper right', labels=['Notnulled', 'Nulled'], fontsize=18, title_fontsize=20)
plt.tight_layout()
plt.savefig(f'results/group_scrubbedVols.png', bbox_inches="tight")
plt.show()
