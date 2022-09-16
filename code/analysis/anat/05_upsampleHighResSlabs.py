'''

Upsampling anatomy and polished segmentation

'''

import os
import glob
import subprocess
# Set data path
DATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-06']

for sub in subs:

    file = glob.glob(f'{DATADIR}/derivatives/{sub}/anat/{sub}_ses-0*_highres-mp2rage_average_uni.nii')[0]
    for i in range(1,4):
        if f'ses-0{i}' in file:
            ses = f'ses-0{i}'

    # Define output folder
    outFolder = f'{DATADIR}/derivatives/{sub}/anat/upsampled'
    # Create folder if not exists
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        print("Subject directory is created")


    # =========================================================================
    # Upsample anatomy

    # Load data
    inFile = f'{DATADIR}/derivatives/{sub}/anat/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected_trunc.nii.gz'
    outFile = f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected_trunc_upsampled.nii.gz'

    # Upsamling
    command = f'c3d '
    command += f'{inFile} '
    command += f'-resample 300x300x300% '
    command += f'-interpolation Cubic '
    command += f'-o {outFile}'

    subprocess.run(command, shell = True)

    # # =========================================================================
    # # Upsample segmentation
    #
    #
    # # Load data
    # inFile = f'{DATADIR}/derivatives/{sub}/anat/sub-06_seg_rim_trunc_polished.nii.gz'
    # outFile = f'{outFolder}/{sub}_seg_rim_trunc_polished_upsampled.nii.gz'
    #
    # # Upsamling
    # command = f'c3d '
    # command += f'{inFile} '
    # command += f'-split -foreach '
    # command += f'-resample 300x300x300% '
    # command += f'-interpolation Cubic '
    # command += f'-smooth 2x2x2vox '
    # command += f'-endfor -merge '
    # command += f'-o {outFile}'
    #
    # subprocess.run(command, shell = True)
