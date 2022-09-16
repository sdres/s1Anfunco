'''

Doing LayNii stuff like:
- layerification
- Defining perimeter
- Flattening

'''

import subprocess
import nibabel as nb
import numpy as np

# Set data path
DATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives'

# Set subjects to work on
subs = ['sub-06']

# Defining control points for all subjects
controlPoints = {'sub-06': [258, 169, 85]}
radii = {'sub-06': 11}

for sub in subs:
    subDir = f'{DATADIR}/{sub}/anat/upsampled'

    # Set rim file
    rimFile = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled.nii.gz'

    # =========================================================================
    # LN2_LAYERS
    command = f'LN2_LAYERS '
    command += f'-rim {rimFile} '
    command += f'-nr_layers 12 '
    command += f'-curvature '
    command += f'-thickness '
    command += f'-streamlines '
    command += '-equivol'

    subprocess.run(command, shell = True)

    # =========================================================================
    # Making control points file

    # Get dimensions of anatomy
    file = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled_midGM_equivol.nii.gz'
    nii = nb.load(file)
    affine = nii.affine
    header = nii.header
    data = nii.get_fdata()

    data[controlPoints[sub][0], controlPoints[sub][1], controlPoints[sub][2]] = 2

    cpFile = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled_midGM_equivol_cp.nii.gz'
    # Save control point
    ni_img = nb.Nifti1Image(cp.astype('int'), affine=affine, header=header)
    nb.save(ni_img, f'{cpFile}')


    # =========================================================================
    # LN2_Multilaterate

    command = f'LN2_MULTILATERATE '
    command += f'-rim {rimFile} '
    command += f'-control_points {cpFile} '
    command += f'-radius {radii[sub]}'

    subprocess.run(command, shell = True)


    # =========================================================================
    # Flattening

    base = f'{subDir}/{sub}_seg_rim_trunc_polished_upsampled'

    command = f'LN2_PATCH_FLATTEN '
    command += f'-coord_uv {base}_UV_coordinates.nii.gz '
    command += f'-coord_d {base}_metric_equivol.nii.gz '
    command += f'-domain {base}_perimeter_chunk.nii.gz '
    command += f'-bins_u 1000 -bins_v 1000 -bins_d 100 '
    command += f'-values {base}_curvature_binned.nii.gz '
    command += f'-voronoi'

    subprocess.run(command, shell = True)
