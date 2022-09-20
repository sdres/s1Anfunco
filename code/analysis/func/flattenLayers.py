import subprocess


subs = ['sub-02', 'sub-05', 'sub-06', 'sub-07', 'sub-09', 'sub-10', 'sub-12', 'sub-15', 'sub-16', 'sub-17', 'sub-18']

for sub in subs:
    rimFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled'
    outFolder = f'/home/sebastian/Desktop/patchFlatten/{sub}_done/upsampled/registeredFunc'


    print(f'flattening layers...')

    inFile = f'{rimFolder}/seg_rim_polished_layers_equivol.nii.gz'
    # print(f'using file: {inFile}')
    command = f'/home/sebastian/git/laynii/LN2_PATCH_FLATTEN '
    command += f'-coord_uv {rimFolder}/seg_rim_polished_UV_coordinates.nii.gz '
    command += f'-coord_d {rimFolder}/seg_rim_polished_metric_equivol.nii.gz '
    command += f'-domain {rimFolder}/seg_rim_polished_perimeter_chunk.nii.gz '
    command += f'-bins_u 1000 -bins_v 1000 -bins_d 100 '
    command += f'-values {inFile} '
    command += f'-voronoi'
    subprocess.run(command, shell=True)


    print('fixing header in layers file...')
    command = f'fslmaths {outFolder}/D2vsAll_BOLD_flat_1000x1000_voronoi.nii '
    command += f'-mul 0 '
    command += f'-add {rimFolder}/seg_rim_polished_layers_equivol_flat_1000x1000_voronoi.nii.gz '
    command += f'{rimFolder}/seg_rim_polished_layers_equivol_flat_1000x1000_voronoi.nii.gz'
    subprocess.run(command, shell=True)
