import nibabel as nb
import numpy as np
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from collections import Counter

subs = ['sub-02', 'sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

# for sub in subs:
#     rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
#
#
#     STREAMLINES = f'{rimFolder}/seg_rim_polished_streamline_vectors.nii.gz'
#     nii2 = nb.load(STREAMLINES)
#     vec_local = np.asarray(nii2.dataobj)
#     vec_local.shape
#
#     vectorImgAnat = np.zeros(vec_local.shape)
#
#     for i in range(3):
#         vectorImgAnat[...,-1]=1
#
#     term1 = np.sqrt(np.sum(vectorImgAnat ** 2., axis=-1)) # [1,0,0]
#     term2 = np.sqrt(np.sum(vec_local ** 2., axis=-1)) # [0,1,0]
#
#     # comparison
#     temp_dot = np.sum(vectorImgAnat * vec_local, axis=-1)
#     temp_angle = np.arccos(temp_dot / (term1 * term2))
#     temp_angle = temp_angle * 180 / np.pi
#     # temp_angle[~idx] = 0
#     temp_angle[np.isnan(temp_angle)] = 0
#     # temp_angle[idx] = 90 - np.abs(temp_angle[idx] - 90)
#     np.percentile(temp_angle, (0,100))
#     # Save
#     outname = os.path.join(rimFolder, "seg_rim_polished_streamline_vectors_localdiff.nii")
#     img = nb.Nifti1Image(temp_angle, affine=nii2.affine, header=nii2.header)
#     nb.save(img, outname)



voxelList = []
angleList = []
subList = []
thicknessList = []
tsnrList = []
kurtList = []
skewList = []
# vasoList = []
# boldList = []

for sub in subs:
    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'

    dataRoot='/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'
    runs = sorted(glob.glob(f'{dataRoot}/{sub}/ses-00*/func/{sub}_ses-00*_task-stim_run-00*_cbv.nii.gz'))
    base = os.path.basename(runs[0]).rsplit('.', 2)[0][:-4]

    # extracting values from perimeter chunk
    mask = f'{rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    maskNii = nb.load(mask)
    maskData = np.asarray(maskNii.dataobj)
    maskIDX = maskData == 1

    anlgeFile = f'{rimFolder}/seg_rim_polished_streamline_vectors_localdiff.nii'
    angleNii = nb.load(anlgeFile)
    angleData = angleNii.get_fdata()
    vecData = angleData[maskIDX]

    thickFile = f'{rimFolder}/seg_rim_polished_thickness.nii.gz'
    thickNii = nb.load(thickFile)
    thickData = thickNii.get_fdata()
    thickness = thickData[maskIDX]

    tsnrFile = f'{rimFolder}/registeredFunc/{base}_VASO_LN_tSNR.nii'
    tsnrNii = nb.load(tsnrFile)
    tsnrData = tsnrNii.get_fdata()
    tsnrMasked = tsnrData[maskIDX]

    kurtFile = f'{rimFolder}/registeredFunc/{base}_VASO_LN_kurt.nii'
    kurtNii = nb.load(kurtFile)
    kurtData = kurtNii.get_fdata()
    kurtMasked = kurtData[maskIDX]

    skewFile = f'{rimFolder}/registeredFunc/{base}_VASO_LN_skew.nii'
    skewNii = nb.load(skewFile)
    skewData = skewNii.get_fdata()
    skewMasked = skewData[maskIDX]

    for i, (val,thick,tsnr, kurt, skew) in enumerate(zip(vecData,thickness,tsnrMasked, kurtMasked, skewMasked)):
        subList.append(sub)
        voxelList.append(i)
        angleList.append(val)
        thicknessList.append(thick)
        tsnrList.append(tsnr)
        kurtList.append(kurt)
        skewList.append(skew)


data = pd.DataFrame({'subject': subList, 'idx': voxelList, 'angle': angleList, 'thickness': thicknessList, 'tSNR': tsnrList, 'kurt': kurtList, 'skew': skewList})


thickCount = Counter()
# thickCount.update(np.round(tmp1['thickness'],decimals=4))
thickCount.update(tmp1['thickness'])
mostCommon = thickCount.most_common()



import math

print((0.1667*math.sqrt(3)*2))

2.4/0.16


digits=['D2','D3','D4']

subList = []
boldList = []
vasoList = []
lenList = []
degDiffList = []
relResList = []
tsnrListBold = []
tsnrListVaso = []

for sub in subs:
    subList.append(sub)

    tmp1 = data.loc[data['subject']==sub]
    medAngle = np.median(tmp1['angle'])

    angles = Counter()
    angles.update(tmp1['angle'].round(decimals=3))
    mostCommon = angles.most_common()[0][0]

    angleDiff = abs(mostCommon-90)
    degDiffList.append(angleDiff)

    medThhickness = np.median(tmp1['thickness'])
    relRes = medThhickness/0.75 # in plane voxel size
    relResList.append(relRes)

    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    dataRoot='/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'

    runs = sorted(glob.glob(f'{dataRoot}/{sub}/ses-00*/func/{sub}_ses-00*_task-stim_run-00*_cbv.nii.gz'))
    base = os.path.basename(runs[0]).rsplit('.', 2)[0][:-4]

    lenList.append(len(runs))

    mask = f'{rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    maskNii = nb.load(mask)
    maskData = np.asarray(maskNii.dataobj)
    maskIDX = maskData == 1

    tsnrFile = f'{rimFolder}/registeredFunc/{base}_VASO_LN_tSNR.nii'
    tsnrNii = nb.load(tsnrFile)
    tsnrData = tsnrNii.get_fdata()
    tsnrMasked = tsnrData[maskIDX]
    tsnrMean = np.mean(tsnrMasked)
    tsnrListVaso.append(tsnrMean)

    tsnrFile = f'{rimFolder}/registeredFunc/{base}_BOLD_intemp_tSNR.nii'
    tsnrNii = nb.load(tsnrFile)
    tsnrData = tsnrNii.get_fdata()
    tsnrMasked = tsnrData[maskIDX]
    tsnrMean = np.mean(tsnrMasked)
    tsnrListBold.append(tsnrMean)

    for modality in ['BOLD','VASO']:
        tmp = []
        for digit in digits:
            vasoFile = f'{rimFolder}/registeredFunc/{digit}vsAll_{modality}.nii'
            vasoNii = nb.load(vasoFile)
            vasoData = vasoNii.get_fdata()
            vasoMasked = vasoData[maskIDX]

            for i, val in enumerate(vasoMasked):
                tmp.append(val)

        tmp = np.array(tmp)
        tmp = np.sort(tmp)[::-1]

        fivePer = round((tmp.shape[0]/100)*5)
        meanAct = np.mean(tmp[:fivePer])

        if modality == 'BOLD':
            boldList.append(meanAct)
        if modality == 'VASO':
            vasoList.append(meanAct)


activatinData = pd.DataFrame({'subject':subList, 'BOLD':boldList, 'VASO': vasoList, 'nrRuns':lenList, 'angDiff':degDiffList, 'relRes': relResList,'motion':motList, 'BOLD tSNR':tsnrListBold, 'VASO tSNR': tsnrListVaso})



stats.pearsonr(activatinData["tSNR"],activatinData['VASO'])

fig, ax = plt.subplots()
sns.scatterplot(data=activatinData, x="angDiff", y="VASO", hue='subject',s=100)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.xlabel('Cortical Thickness (mm)', fontsize= 24)
plt.ylabel('5% highest VASO (z)', fontsize= 24)
plt.xlabel('Most Common Angular Error (deg.)', fontsize= 24)
# plt.xticks(np.arange(17,24),fontsize=20)
plt.yticks(fontsize=20)
# plt.yticks([],fontsize=20)
plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/VASO-angDiff_jmostCommon.png',bbox_inches='tight')
plt.show()



fig, ax = plt.subplots()
sns.scatterplot(data=activatinData, x="BOLD", y="VASO", hue='subject',s=100)
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# plt.xlabel('Cortical Thickness (mm)', fontsize= 24)
plt.ylabel('mean of\n5% highest VASO (z)', fontsize= 24)
plt.xlabel('mean of\n5% highest BOLD (z)', fontsize= 24)
# plt.xticks(np.arange(17,24),fontsize=20)
plt.yticks(fontsize=20)
# plt.yticks([],fontsize=20)
plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/VASO-BOLD.png',bbox_inches='tight')
plt.show()


for modality in ['BOLD', 'VASO']:
    fig, ax = plt.subplots()
    sns.scatterplot(data=activatinData, x=f"{modality} tSNR", y=modality, hue='subject',s=100)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.ylabel(f'mean of\n5% highest zscore', fontsize= 24)
    plt.xlabel('tSNR', fontsize= 24)
    # plt.xticks(np.arange(17,24),fontsize=20)
    plt.yticks(fontsize=20)
    # plt.yticks([],fontsize=20)
    plt.title(modality, fontsize=24)
    plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/{modality}-tSNR.png',bbox_inches='tight')
    plt.show()




fig, ax = plt.subplots(1,2)
sns.boxplot(ax= ax[0], data = activatinData, x='nrRuns',y='BOLD')
sns.boxplot(ax= ax[1],data = activatinData, x='nrRuns',y='VASO')
ax[0].set_title('BOLD', fontsize=18)
ax[1].set_title('VASO', fontsize=18)
ax[0].set_ylabel('mean of\n5% highest z-score', fontsize= 24)
ax[1].set_ylabel('', fontsize= 18)
ax[0].set_xlabel('Nr. Runs', fontsize= 18)
ax[1].set_xlabel('Nr. Runs', fontsize= 18)
plt.savefig('/home/sebastian/Desktop/patchFlatten/screenshots/vasoBOLDNrRuns.png',bbox_inches='tight')
plt.show()



fig = plt.figure()
sns.kdeplot(data=activatinData, x='relRes', hue='subject',legend=True)
# med = np.median(tmp['angle'])
# plt.axvline(med, color='black',label=f'median: {med}')
# plt.title(sub)
# plt.legend()
# plt.xaxis('degrees')
plt.show()


from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter


fig, ax = plt.subplots()
sns.kdeplot(data=data,x='thickness',legend=True, linewidth=2,hue='subject')
# sns.histplot(data=tmp,x='thickness')
med = np.median(data['thickness'])

# angles = Counter()
# angles.update(tmp['angle'].round(decimals=3))
# mostCommon = angles.most_common()[0][0]

# plt.axvline(mostCommon, color='red',label=f'most common: {mostCommon}')

plt.axvline(med, color='black',label=f'median: {med:1f}')
ax.xaxis.set_major_locator(MaxNLocator(integer=False))

# plt.title(sub, fontsize= 24)
plt.legend()
plt.xlabel('Cortical Thickness (mm)', fontsize= 24)
plt.ylabel('Voxel Count', fontsize= 24)
plt.xticks(np.arange(1,4.5,0.5),fontsize=20)
plt.yticks([],fontsize=20)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/allSubsIndivid_cortThickness',bbox_inches='tight')
plt.show()




for sub in subs:
    tmp = data.loc[data['subject']==sub]

    fig, ax = plt.subplots()
    sns.kdeplot(data=tmp,x='thickness',legend=True, linewidth=2)
    # sns.histplot(data=tmp,x='thickness')
    med = np.median(tmp['thickness'])

    # angles = Counter()
    # angles.update(tmp['angle'].round(decimals=3))
    # mostCommon = angles.most_common()[0][0]

    # plt.axvline(mostCommon, color='red',label=f'most common: {mostCommon}')

    plt.axvline(med, color='black',label=f'median: {med:1f}')
    ax.xaxis.set_major_locator(MaxNLocator(integer=False))

    plt.title(sub, fontsize= 24)
    plt.legend()
    plt.xlabel('Cortical Thickness (mm)', fontsize= 24)
    plt.ylabel('Voxel Count', fontsize= 24)
    plt.xticks(np.arange(1,4.5,0.5),fontsize=20)
    plt.yticks([],fontsize=20)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/{sub}_cortThickness',bbox_inches='tight')
    plt.show()


for sub in subs:
    tmp = data.loc[data['subject']==sub]

    fig, ax = plt.subplots()
    ax = sns.kdeplot(data=tmp,x='angle',legend=True, linewidth=2)
    med = np.median(tmp['angle'])

    # y = ax.lines[0].get_ydata() # Get the y data of the distribution
    # maxid = np.max(y)

    angles = Counter()
    angles.update(tmp['angle'].round(decimals=3))
    mostCommon = angles.most_common()[0][0]

    plt.axvline(mostCommon, color='red',label=f'most common: {mostCommon}')

    plt.axvline(med, color='black',label=f'median: {med:1f}')
    ax.xaxis.set_major_locator(MaxNLocator(integer=False))

    plt.title(sub,fontsize=24)
    plt.legend()
    plt.xlabel('LocalGM-Slab Angle (deg)', fontsize= 24)
    plt.ylabel('Voxel Count', fontsize= 24)
    plt.xticks(np.arange(10,141,10),fontsize=10)
    plt.yticks([],fontsize=20)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))

    plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/{sub}_angleDiff.png',bbox_inches='tight')
    plt.show()

### Get motion



motList = []
for sub in subs:
    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    dataRoot='/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'
    runs = sorted(glob.glob(f'{dataRoot}/{sub}/ses-00*/func/{sub}_ses-00*_task-stim_run-00*_cbv.nii.gz'))
    # base = os.path.basename(runs[0]).rsplit('.', 2)[0][:-4]
    tmpSub = []
    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        motionData = pd.read_csv(f'{dataRoot}/derivatives/{sub}/func/motionParameters/{base}_FDs.csv')
        for modality in ['BOLD','VASO']:
            tmp = motionData.loc[(motionData['modality']==modality)&(motionData['run']==base)]
            tmpSub.append(np.sum(tmp['FD']))
    tmpSub = np.sum(tmpSub)
    tmpSub = tmpSub/2

    tmpSub = tmpSub/len(runs)
    motList.append(tmpSub)








digits=['D2','D3','D4']

subList = []
modalityList = []
actList = []
digitList = []

for sub in subs:


    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    dataRoot='/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS'

    mask = f'{rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    maskNii = nb.load(mask)
    maskData = np.asarray(maskNii.dataobj)
    maskIDX = maskData == 1

    for modality in ['BOLD','VASO']:

        for digit in digits:
            tmp = []
            vasoFile = f'{rimFolder}/registeredFunc/{digit}vsAll_{modality}.nii'
            vasoNii = nb.load(vasoFile)
            vasoData = vasoNii.get_fdata()
            vasoMasked = vasoData[maskIDX]

            for i, val in enumerate(vasoMasked):
                tmp.append(val)

            tmp = np.array(tmp)
            tmp = np.sort(tmp)[::-1]

            fivePer = round((tmp.shape[0]/100)*5)
            meanAct = np.mean(tmp[:fivePer])

            modalityList.append(modality)
            actList.append(meanAct)
            digitList.append(digit)
            subList.append(sub)



activatinDataPerDigit = pd.DataFrame({'subject':subList, 'modality':modalityList, 'activation': actList, 'digit':digitList})

for modality in ['BOLD','VASO']:
    sns.barplot(data=activatinDataPerDigit.loc[activatinDataPerDigit['modality']==modality],x='subject',y='activation', hue='digit')
    plt.title(modality, fontsize=24)
    plt.xlabel('')
    plt.ylabel('mean of\n5% highest z-score', fontsize=18)
    plt.savefig(f'/home/sebastian/Desktop/patchFlatten/screenshots/digitActivations_{modality}.png',bbox_inches='tight')
    plt.show()



#
