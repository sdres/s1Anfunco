'''

Running first level GLM in FSL using Nilearn

'''

import nibabel as nb
import nilearn
import numpy as np
from nilearn.glm.first_level import make_first_level_design_matrix
import glob
import pandas as pd
from nilearn.glm.first_level import FirstLevelModel


ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'


tr = 1.9295

subs = ['sub-06']

drift_model = 'Cosine'  # We use a discrete cosine transform to model signal drifts.
high_pass = .01  # The cutoff for the drift model is 0.01 Hz.
hrf_model = 'spm'  # The hemodynamic response function is the SPM canonical one.


for sub in subs:

    funcDir = f'{ROOT}/derivatives/{sub}/func'
    # make folder to dump statistocal maps
    statFolder = f'{funcDir}/statMaps'
    if not os.path.exists(statFolder):
        os.makedirs(statFolder)
        print("Statmap directory is created")

    for ses in ['ses-01']:

        runs = sorted(glob.glob(f'{ROOT}/{sub}/{ses}/func/{sub}_{ses}_task-stim*run-00*_cbv.nii.gz'))


        for run in runs:
            base = os.path.basename(run).rsplit('.', 2)[0][:-4]

            for modality in ['VASO', 'BOLD']:

                niiFile = f'{funcDir}/{base}_{modality}.nii.gz'
                nii = nb.load(niiFile)
                data = nii.get_fdata()
                nVols = data.shape[-1]
                frame_times = np.arange(nVols) * tr

                events = pd.read_csv(f'{ROOT}/{sub}/{ses}/func/{base}_events.tsv', sep = ' ')

                design_matrix = make_first_level_design_matrix(
                    frame_times,
                    events,
                    hrf_model=hrf_model,
                    drift_model = None,
                    high_pass= high_pass
                    )

                contrast_matrix = np.eye(design_matrix.shape[1])
                basic_contrasts = dict([(column, contrast_matrix[i])
                            for i, column in enumerate(design_matrix.columns)])

                if modality == 'BOLD':
                    contrasts = {
                        'D2VsAll': + basic_contrasts['D2'] - basic_contrasts['D3']/2 - basic_contrasts['D4']/2,
                        'D3VsAll': - basic_contrasts['D2']/2 + basic_contrasts['D3'] - basic_contrasts['D4']/2,
                        'D4VsAll': - basic_contrasts['D2']/2 - basic_contrasts['D3']/2 + basic_contrasts['D4'],
                        'D2VsRest': + basic_contrasts['D2'],
                        'D3VsRest': + basic_contrasts['D3'],
                        'D4VsRest': + basic_contrasts['D4']
                        }

                if modality == 'VASO':
                    contrasts = {
                        'D2VsAll': - basic_contrasts['D2'] + basic_contrasts['D3']/2 + basic_contrasts['D4']/2,
                        'D3VsAll': + basic_contrasts['D2']/2 - basic_contrasts['D3'] + basic_contrasts['D4']/2,
                        'D4VsAll': + basic_contrasts['D2']/2 + basic_contrasts['D3']/2 - basic_contrasts['D4'],
                        'D2VsRest': - basic_contrasts['D2'],
                        'D3VsRest': - basic_contrasts['D3'],
                        'D4VsRest': - basic_contrasts['D4']
                        }

                fmri_glm = FirstLevelModel(mask_img = False, drift_model=None)
                fmri_glm = fmri_glm.fit(nii, design_matrices = design_matrix)

                # Iterate on contrasts
                for contrast_id, contrast_val in contrasts.items():
                    # compute the contrasts
                    z_map = fmri_glm.compute_contrast(
                        contrast_val, output_type='z_score')
                    nb.save(z_map, f'{statFolder}/{base}_{modality}_{contrast_id}.nii')
