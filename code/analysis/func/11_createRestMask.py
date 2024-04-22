"""Create mask to minimaze storage demands of resting state analysis"""

import subprocess
import os
import nibabel as nb
import numpy as np
from scipy.ndimage import morphology, generate_binary_structure
from nipype.interfaces.fsl import ExtractROI


# Set rough outlines of the post-central sulcus (visually defined)
boundariesDict = {'sub-14': [80, 170, 75, 115, 0, 60],
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
ROOT = f'/Users/sebastiandresbach/data/s1Anfunco/Nifti'
subs = ['sub-14', 'sub-15', 'sub-16', 'sub-17', 'sub-18']
for sub in subs:

    # ===================================
    # Get brain mask from brain extraction
    inFolder = f'/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives_old/manualSteps/{sub}'
    outFolder = f'/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives/{sub}/anat'
    inFile = f'{sub}_highres-mp2rage_uni_average_N4corrected_brainMask_opening_closing.nii.gz'
    outFile = f'{sub}_brainMask.nii.gz'
    command = f'cp {inFolder}/{inFile} {outFolder}/{outFile}'
    subprocess.run(command, shell=True)

    # ===================================
    # Shrink brain mask to make sure that only brain and CSF is included

    nii = nb.load(f"{outFolder}/{outFile}")
    data = np.asarray(nii.dataobj)
    data = np.pad(data, 1, mode="reflect")  # to prevent data edge artifacts

    struct = generate_binary_structure(3, 1)  # 1 jump neighbourbhood

    data = morphology.binary_erosion(data, structure=struct, iterations=10)

    SUFFIX = "erode"
    basename, ext = nii.get_filename().split(os.extsep, 1)
    out = nb.Nifti1Image(data.astype(int), header=nii.header, affine=nii.affine)
    nb.save(out, "{}_{}.{}".format(basename, SUFFIX, ext))

    # ===================================
    # Crop mask
    inFile = "{}_{}.{}".format(basename, SUFFIX, ext)
    outFile = "{}_{}_trunc.{}".format(basename, SUFFIX, ext)

    # Apparently, the fslroi wrapper in nipype wants an existing file as output filename (?).
    # Therefore, I will create one with the same dimensions.
    tmpFile = nb.load(inFile)

    tmp = nb.Nifti1Image(tmpFile.get_fdata(), header=tmpFile.header, affine=tmpFile.affine)

    nb.save(tmp, outFile)

    fslroi = ExtractROI(in_file=inFile,
                        roi_file=outFile,
                        x_min=boundariesDict[sub][0], x_size=boundariesDict[sub][1],
                        y_min=boundariesDict[sub][2], y_size=boundariesDict[sub][3],
                        z_min=boundariesDict[sub][4], z_size=boundariesDict[sub][5]
                        )

    out = fslroi.run()

    # ===================================
    # Resample

    command = f'c3d '
    command += "{}_{}_trunc.{} ".format(basename, SUFFIX, ext)
    command += f'-resample 300x300x300% '
    command += f'-interpolation NearestNeighbor '
    command += "-o {}_{}_trunc_scaled.{} ".format(basename, SUFFIX, ext)

    subprocess.run(command, shell=True)

    command = f'fslmaths {basename}_{SUFFIX}_trunc_scaled.{ext} -bin {basename}_{SUFFIX}_trunc_scaled_bin.{ext}'
    subprocess.run(command, shell=True)

    # ===================================
    # intersect with CSF

    segFolder = f'/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives/{sub}/anat/upsampled'
    segFile = f'{segFolder}/seg_rim_polished.nii.gz'
    maskFile = f'{basename}_erode_trunc_scaled_bin.nii.gz'

    maskNii = nb.load(maskFile)
    maskData = maskNii.get_fdata()
    segNii = nb.load(segFile)
    segData = segNii.get_fdata()

    csf = segData == 1

    intersect = csf * maskData

    new = nb.Nifti1Image(intersect.astype(int), header=maskNii.header, affine=maskNii.affine)
    nb.save(new, f'/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives/{sub}/anat/{sub}_csf.nii.gz')

    anatFolder = f'{ROOT}/derivatives/{sub}/anat/upsampled'

    # Add perimeter
    periFile = f'{anatFolder}/seg_rim_polished_perimeter_chunk.nii.gz'
    periNii = nb.load(periFile)
    periData = periNii.get_fdata()
    perihead = periNii.header
    periAff = periNii.affine

    perNoBorderData = np.where(periData == 1, 1, 0)

    csfPeri = intersect + perNoBorderData

    img = nb.Nifti1Image(csfPeri.astype('int'), header=perihead, affine=periAff)
    nb.save(img, f'{anatFolder}/{sub}_csfPlusPeri.nii.gz')
