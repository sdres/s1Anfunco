"""

Doing LayNii stuff like:
- layerification
- Defining perimeter
- Flattening

"""

import subprocess
from skimage.measure import label
import nibabel as nb
import numpy as np

# Set data path
DATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'

# Set subjects to work on
subs = ['sub-12']
# Set sessions to work on
sessions = ['ses-01']

# Defining control points for all subjects
controlPoints = {'sub-05': [281, 209, 120],
                 'sub-06': [263, 170, 80],
                 'sub-07': [275, 206, 100],
                 'sub-10': [335, 174, 105],
                 'sub-12': [247, 175, 110],
                 'sub-15': [241, 178, 105],
                 'sub-16': [259, 190, 135],
                 'sub-17': [212, 207, 90],
                 'sub-18': [258, 159, 100]}  # 261, 183, 90, 258, 180, 90

radii = {'sub-05': 11,
         'sub-06': 12,
         'sub-07': 11,
         'sub-09': 11,
         'sub-10': 11,
         'sub-12': 11,
         'sub-15': 11,
         'sub-16': 11,
         'sub-17': 11,
         'sub-18': 11}


for sub in subs:
    subDir = f'{DATADIR}/{sub}/anat/upsampled'

    # Set rim file
    # rimFile = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled.nii.gz'
    rimFile = f'{subDir}/seg_rim_polished.nii.gz'

    # =========================================================================
    # LN2_LAYERS
    command = f'LN2_LAYERS '
    command += f'-rim {rimFile} '
    command += f'-nr_layers 11 '
    command += f'-curvature '
    command += f'-thickness '
    command += f'-streamlines '
    command += '-equivol'

    subprocess.run(command, shell=True)

    # =========================================================================
    # Making control points file

    # Get dimensions of anatomy
    file = f'{subDir}/seg_rim_polished_midGM_equivol.nii.gz'
    nii = nb.load(file)
    affine = nii.affine
    header = nii.header
    cp = nii.get_fdata()

    cp[controlPoints[sub][0], controlPoints[sub][1], controlPoints[sub][2]] = 2

    cpFile = f'{subDir}/seg_rim_polished_midGM_equivol_cp.nii.gz'
    # Save control point
    ni_img = nb.Nifti1Image(cp.astype('int'), affine=affine, header=header)
    nb.save(ni_img, f'{cpFile}')

    # =========================================================================
    # LN2_MULTILATERATE

    command = f'LN2_MULTILATERATE '
    command += f'-rim {rimFile} '
    command += f'-control_points {cpFile} '
    command += f'-radius {radii[sub]}'

    subprocess.run(command, shell=True)

    # Remove floating bits
    # periFile = f'{subDir}/seg_rim_polished_perimeter_chunk.nii.gz'
    # mapNii = nb.load(periFile)
    # periData = mapNii.get_fdata()
    # data = label(periData, connectivity=1)
    # labels, counts = np.unique(data, return_counts=True)
    # largestCluster = np.argmax(counts[1:]) + 1
    #
    # tmp = data == largestCluster
    # img = nb.Nifti1Image(tmp, affine=mapNii.affine, header=mapNii.header)
    # nb.save(img, f'{subDir}/seg_rim_polished_perimeter_chunk.nii.gz')

    # =========================================================================
    # Flattening

    # base = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled'
    #
    # command = f'LN2_PATCH_FLATTEN '
    # command += f'-coord_uv {base}_UV_coordinates.nii.gz '
    # command += f'-coord_d {base}_metric_equivol.nii.gz '
    # command += f'-domain {base}_perimeter_chunk.nii.gz '
    # command += f'-bins_u 1000 -bins_v 1000 -bins_d 100 '
    # command += f'-values {base}_curvature_binned.nii.gz '
    # command += f'-voronoi'
    #
    # subprocess.run(command, shell = True)
