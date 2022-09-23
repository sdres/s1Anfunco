'''

Plotting motion traces for visual inspection.

'''




import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import numpy as np

subs = ['sub-06']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'



plt.style.use('dark_background')


for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))


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

                for i, val in enumerate(tmp, start = 1):
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

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))


    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]

        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'


        data = pd.read_csv(os.path.join(motionDir, f'{base}_motionSummary.csv'))

        fig, axes = plt.subplots(1, 2,sharex=True,figsize=(30,6))
        plt.suptitle(f'{base} Motion Summary', fontsize=24)


        for i, type in enumerate(['T', 'R']):

            if type == 'T':
                legend = False
            if type == 'T':
                legend = True

            sns.lineplot(ax=axes[i], x='volume',y='Motion', data=data.loc[(data['name'].str.contains(type)) & (data['name'].str.contains('nulled'))], hue='name', palette = 'Set1',linewidth = width,legend=legend)
            sns.lineplot(ax=axes[i], x='volume',y='Motion', data=data.loc[(data['name'].str.contains(type)) & (data['name'].str.contains('notnulled'))], hue='name', palette = 'Set2',linewidth = width,legend=legend)


        # sns.lineplot(ax=axes[1], x='volume',y='Motion',data=data.loc[(data['name'].str.contains('R')) & (data['name'].str.contains('nulled'))],hue='name', palette = 'Set1',linewidth = width)
        # sns.lineplot(ax=axes[1], x='volume',y='Motion',data=data.loc[(data['name'].str.contains('R')) & (data['name'].str.contains('notnulled'))],hue='name', palette = 'Set2',linewidth = width)


        axes[0].set_ylabel("Translation [mm]", fontsize=24)
        axes[1].set_ylabel("Rotation [radians]", fontsize=24)

        # Colors are the same for both plots so we only need one legend
        axes[1].legend(fontsize=20, loc='center left', bbox_to_anchor=(1, 0.5))

        for j in range(2):
            axes[j].tick_params(axis='both', labelsize=20)
            axes[j].set_xlabel("Volume", fontsize=24)

        plt.savefig(f'{outDir}/{base}_motion.jpg', bbox_inches = 'tight', pad_inches = 0)
        plt.show()




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


        FDs = pd.read_csv(os.path.join(motionDir, f'{base}_FDs.csv'))

        fig = plt.figure(figsize=(20,5))
        sns.lineplot(data=FDs, x='volume', y='FD',hue='modality',linewidth = width, palette=v1Palette)

        if np.max(FDs['FD']) < 0.88:
            plt.ylim(0,1)


        plt.axhline(0.88, color='gray', linestyle='--')
        plt.ylabel('FD [mm]', fontsize=24)
        plt.xlabel('Volume', fontsize=24)
        plt.title(f"{base}", fontsize=24, pad=20)
        plt.legend(fontsize=20)


        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.savefig(f'{outDir}/{base}_FDs.png', bbox_inches = 'tight', pad_inches = 0)
        plt.show()
