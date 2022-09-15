'''

Cropping anatomical images to decrease number of voxels (and therefore
processing demands)

'''
from nipype.interfaces.fsl import ExtractROI
import nibabel as nb
import numpy as np

# Set data path
DATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-05']
# Set sessions to work on
sessions = ['ses-02']

# Set rough outlines of the post-central sulcus (visually defined)
boundariesDict = {'sub-02': [80, 170, 75, 115, 0, 60],
                  'sub-05': [110, 160, 60, 125, 0, 60],
                  'sub-06': [100, 160, 70, 90, 0, 60],
                  'sub-07': [100, 180, 70, 120, 0, 60],
                  'sub-09': [75, 175, 110, 90, 0, 60],
                  'sub-10': [90, 175, 80, 100, 0, 60],
                  'sub-12': [90, 175, 60, 110, 0, 60],
                  'sub-15': [85, 165, 80, 110, 0, 60],
                  'sub-16': [90, 160, 65, 125, 0, 60],
                  'sub-17': [80, 155, 95, 115, 0, 60],
                  'sub-18': [65, 170, 90, 115, 0, 60]
                 }


for sub in subs:
    for ses in sessions:
        # Define output folder
        outFolder = f'{DATADIR}/derivatives/{sub}/anat'


        inFile = f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected.nii'
        outFile = f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected_trunc.nii.gz'

        # Apparently, the fslroi wrapper in nipype wants an existing file as output filename (?). Therefore, I will create one with the same dimensions.
        tmpFile = nb.load(inFile)

        tmp = nb.Nifti1Image(tmpFile.get_fdata(), header=tmpFile.header, affine=tmpFile.affine)

        nb.save(tmp, outFile)

        fslroi = ExtractROI(in_file = inFile,
                            roi_file = outFile,
                            x_min = boundariesDict[sub][0], x_size = boundariesDict[sub][1],
                            y_min = boundariesDict[sub][2], y_size = boundariesDict[sub][3],
                            z_min = boundariesDict[sub][4], z_size = boundariesDict[sub][5]
                           )

        out = fslroi.run()
