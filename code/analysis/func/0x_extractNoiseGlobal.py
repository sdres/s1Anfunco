import nibabel as nb
import glob
import numpy as np
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

root = '/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'

subs = sorted(glob.glob(f'{root}/derivatives/sub*'))
subs = [fileName.split('/')[-1] for fileName in subs]
# test whether all subnames are there
# subs

# initialte lists that will contian the noise measures
subData = []
measureData = []
modalityData = []
runData = []
dataData = []

modalities = ['BOLD_intemp_noTrunc', 'noTrunc_VASO_LN']
measures = ['tSNR', 'skew', 'kurt']

for sub in subs:

    runs = sorted(glob.glob(f'{root}/{sub}/ses-00*/func/{sub}_ses-00*_task-*stim*run-00*_cbv.nii.gz'))

    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        mask = nb.load(f'{root}/derivatives/{sub}/func/{sub}_noiseMask.nii').get_fdata()
        idx = mask == 1

        for measure in measures:
            for modality in modalities:
                data = nb.load(f'{root}/derivatives/{sub}/func/{base}_{modality}_{measure}.nii').get_fdata()[:,:,:,0]
                data = data[idx]
                voxCount = data.size

                if 'BOLD' in modality:
                    modality = 'BOLD'
                if 'VASO' in modality:
                    modality = 'VASO'
                for item in data:

                    subData.append(sub)
                    measureData.append(measure)
                    modalityData.append(modality[:4])
                    runData.append(base)
                    dataData.append(item)

    df = pd.DataFrame({'subject':subData, 'measure':measureData, 'modality':modalityData, 'run':runData, 'data':dataData})




    # g = sns.FacetGrid(df.loc[df['subject']==sub], row="measure", col='modality',sharex=False, sharey=False)
    # g.map(sns.kdeplot, "data")
    # plt.title(sub)
    # plt.show()



voxCount
for sub in subs:
    subtmp = df.loc[(df['modality']=='VASO')&(df['subject']==sub)]
    fig, axes = plt.subplots(1,3,sharey=True)

    for i, measure in enumerate(measures):
        sns.kdeplot(data=subtmp.loc[subtmp['measure']==measure], x='data', hue='run',ax=axes[i])
        # place the legend outside the figure/plot
        axes[i].legend(bbox_to_anchor=(1.1, 1.05),fontsize=20)
        axes[i].set_title(measure, fontsize=20)
    fig.suptitle(sub, fontsize=30)
    fig.tight_layout()
    plt.show()

for measure in measures:
    fig, ax = plt.subplots()
    ax = sns.barplot(data=df.loc[(df['modality']=='VASO')&(df['measure']==measure)],y='data', x='subject')
    plt.title(measure, fontsize=20)
    plt.xticks(rotation=90)
    plt.show()

|import pingouin as pg

aov = pg.anova(dv='data', between='subject', data=df.loc[(df['modality']=='VASO')&(df['measure']=='tSNR')])
aov
pg.pairwise_ttests(dv='data', between='subject', data=df.loc[(df['modality']=='VASO')&(df['measure']=='tSNR')]).round(3)
