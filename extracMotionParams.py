# basic
import glob
import itertools
import pandas as pd
import numpy as np
import os

# plotting
from nilearn import plotting
from scipy import ndimage
import matplotlib.cm as cm
import matplotlib.image as mpimg
import seaborn as sns
import matplotlib.pyplot as plt

# neuroimaging
from nilearn.connectome import ConnectivityMeasure
import nibabel as nb
from nibabel import load, save, Nifti1Image
import nipype.interfaces.fsl as fsl
import nilearn
import ants

# misc
from scipy.ndimage import morphology
from scipy import stats
from scipy import interpolate
import nipype
import logging
logging.getLogger('nipype.interfaces').setLevel(0)

root = '/media/sebastian/Data/S1ANFUNCO/raw_data'

stimRuns = sorted(glob.glob('/media/sebastian/Data/S1ANFUNCO/raw_data/S1ANFUNCO_BIDS/sub-*/ses-00*/func/sub-*_ses-00*_task-stim*_cbv.nii.gz'))
runNumberList = [i[-18:-11] for i in stimRuns]
subList = [i[-43:-37] for i in stimRuns]
outFolders =  sorted(glob.glob('/media/sebastian/Data/S1ANFUNCO/raw_data/derivatives/sub-*/func/stim/motionParameters/run-00*'))

lstsubMot  = []
lstsubMot_Nme  = []
lstTR_sub = []
sub_list =[]
modalityList = []
runsList = []

modalities = ["n", 'nn']
modalitiesDict = {'n': 'VASO', 'nn': 'BOLD'}


for sub in subList:
    for run in runNumberList:

        for modality in modalities:

            print(sub)
            affinemats = sorted(glob.glob(root + f'/derivatives/{sub}/func/stim/motionParameters/{run}/{modality}_motion/vol_*_0GenericAffine_FSL.mat'))
            print(len(affinemats))

            volids=[]
            for n in range(len(affinemats)):
                volids.append(n+1000)

            for volid in volids:

                # Current mats
                currMats = root + f'/derivatives/{sub}/func/stim/motionParameters/{run}/{modality}_motion/vol_{volid}_0GenericAffine_FSL.mat'

                tmp = fsl.AvScale(all_param=True,mat_file=currMats);

                tmpReadout = tmp.run();

                # Get the rotations (in rads) and translations (in mm) per volume

                aryTmpMot = list(itertools.chain.from_iterable([tmpReadout.outputs.translations,tmpReadout.outputs.rot_angles]));




                # Save the roation and translations
                lstsubMot.append(aryTmpMot)
                lstTR_sub.append([int(volid)+1-1000 for i in range(6)])
                lstsubMot_Nme.append([f'TX {modalitiesDict[modality]}',f'TY {modalitiesDict[modality]}',f'TZ {modalitiesDict[modality]}',f'RX {modalitiesDict[modality]}',f'RY {modalitiesDict[modality]}',f'RZ {modalitiesDict[modality]}'])
                sub_list.append([str(sub) for i in range(6)])
                modalityList.append([str(modalitiesDict[modality]) for i in range(6)])
                runsList.append([str(run) for i in range(6)])




aryCurr = np.array(lstsubMot)
aryCurr_Ses =  aryCurr.reshape((aryCurr.size,-1))
aryCurr_TR = np.array(lstTR_sub)
aryCurr_TR_Ses = aryCurr_TR.reshape((aryCurr_TR.size,-1))
aryCurr_Nme = np.array(lstsubMot_Nme)
aryCurr_Nme_Ses = aryCurr_Nme.reshape((aryCurr_Nme.size,-1))
aryIdx = np.arange(1,len(aryCurr_Nme_Ses)+1)

aryCurr_mod = np.array(modalityList)
aryCurr_mod = aryCurr_mod.reshape((aryCurr_mod.size,-1))

aryCurr_subname = np.array(sub_list)
aryCurr_subs =  aryCurr_subname.reshape((aryCurr_subname.size,-1))

aryCurr_runsList = np.array(runsList)
aryCurr_runs =  aryCurr_runsList.reshape((aryCurr_runsList.size,-1))



data_dict = {
    'subject': aryCurr_subs[:,0],
    'Time/TR': aryCurr_TR_Ses[:,0],
    'Motion_Name': aryCurr_Nme_Ses[:,0],
    'Motion': aryCurr_Ses[:,0],
    'idx':aryIdx,
    'modality': aryCurr_mod[:,0],
    'run': aryCurr_runs[:,0]}

pd_ses = pd.DataFrame(data=data_dict)
pd_ses.to_csv(root + '/derivatives/quality_assessmen/motionParametersStim.csv', index=False)
