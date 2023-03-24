"""

Cropping functional images to decrease number of voxels (and therefore
processing demands)

"""

from nipype.interfaces.fsl import ExtractROI
import glob
import os

# Set data path
ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-12']

# Define boundaries for each participant
boundariesDict = {'sub-02': [30, 110, 50, 70, 0, 21],
                  'sub-05': [10, 95, 40, 85, 0, 21],
                  'sub-06': [10, 110, 50, 70, 0, 21],
                  'sub-07': [5, 120, 30, 95, 0, 21],
                  'sub-09': [25, 105, 80, 60, 0, 21],
                  'sub-10': [10, 110, 40, 90, 0, 21],
                  'sub-12': [20, 110, 40, 90, 0, 21],
                  'sub-15': [20, 115, 40, 90, 0, 21],
                  'sub-16': [20, 110, 35, 95, 0, 21],
                  'sub-17': [30, 105, 70, 70, 0, 21],
                  'sub-18': [30, 115, 40, 90, 0, 21]
                  }

# Loop over participants
for sub in subs:

    # Look for individual runs (containing both nulled and notnulled images)
    runs = sorted(glob.glob(f'{ROOT}/{sub}/ses-0*/func/{sub}_ses-0*_task-*run-00*_cbv.nii.gz'))

    # Set session
    for i in range(1, 4):
        if f'ses-0{i}' in runs[0]:
            ses = f'ses-0{i}'

    # Loop over runs
    for run in runs:

        # Get basename of run
        base = os.path.basename(run).rsplit('.', 2)[0][:-4]
        print(f'Processing run {base}')

        # Define output folder
        outFolder = f'{ROOT}/derivatives/{sub}/func'

        # Loop over modalities
        for modality in ['nulled', 'notnulled']:

            inFile = f'{outFolder}/{base}_{modality}_moco.nii.gz'
            outFile = f'{outFolder}/{base}_{modality}_moco_trunc.nii.gz'

            # Prepare command
            fslroi = ExtractROI(in_file=inFile,
                                roi_file=outFile,
                                x_min=boundariesDict[sub][0], x_size=boundariesDict[sub][1],
                                y_min=boundariesDict[sub][2], y_size=boundariesDict[sub][3],
                                z_min=boundariesDict[sub][4], z_size=boundariesDict[sub][5]
                                )
            # Run command
            out = fslroi.run()

        # Do the same for t1w images
        inFile = f'{outFolder}/{sub}_{ses}_T1w.nii'
        outFile = f'{outFolder}/{sub}_{ses}_T1w_trunc.nii.gz'

        fslroi = ExtractROI(in_file=inFile,
                            roi_file=outFile,
                            x_min=boundariesDict[sub][0], x_size=boundariesDict[sub][1],
                            y_min=boundariesDict[sub][2], y_size=boundariesDict[sub][3],
                            z_min=boundariesDict[sub][4], z_size=boundariesDict[sub][5]
                            )

        out = fslroi.run()
