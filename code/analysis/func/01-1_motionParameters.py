'''

Read and plot motion traces

'''

import ants
import os
import glob
import nibabel as nb
import numpy as np
import subprocess
from IPython.display import clear_output
import nipype.interfaces.fsl as fsl
import itertools
import pandas as pd


def my_ants_affine_to_distance(affine, unit):

    dx, dy, dz = affine[9:]

    rot_x = np.arcsin(affine[6])
    cos_rot_x = np.cos(rot_x)
    rot_y = np.arctan2(affine[7] / cos_rot_x, affine[8] / cos_rot_x)
    rot_z = np.arctan2(affine[3] / cos_rot_x, affine[0] / cos_rot_x)

    if unit == 'deg':
        deg = np.degrees
        R = np.array([deg(rot_x), deg(rot_y), deg(rot_z)])
    if unit == 'rad':
        R = np.array([rot_x, rot_y, rot_z])

    T = np.array([dx, dy, dz])

    return T, R


subs = ['sub-06']
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'

    # look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))


    for run in runs:
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')
        # Set folder where motion traces were dumped
        motionDir = f'{funcDir}/motionParameters/{base}'

        for modality in ['nulled', 'notnulled']:
            mats = sorted(glob.glob(f'{motionDir}/{base}_{modality}_vol*'))

            Tr = []; Rt = []

            for i, mat in enumerate(mats):

                localtxp = ants.read_transform(mat)
                affine = localtxp.parameters

                T, R = my_ants_affine_to_distance(affine, 'rad')

                Tr.append(T)
                Rt.append(R)


            # // Save motion traces intra-run as .csv
            Tr = np.asarray(Tr)
            Rt = np.asarray(Rt)

            data_dict = {
            'Tx': Tr[:, 0],
            'Ty': Tr[:, 1],
            'Tz': Tr[:, 2],
            'Rx': Rt[:, 0],
            'Ry': Rt[:, 1],
            'Rz': Rt[:, 2]
            }

            pd_ses = pd.DataFrame(data=data_dict)
            pd_ses.to_csv(os.path.join(motionDir, f'{base}_{modality}_motionRegressors.csv'), index=False)