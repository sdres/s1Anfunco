'''

Registering high resolution anatomical images to boost SNR

'''

import glob
import ants
import numpy as np
import os
import nibabel as nb


# Set data path
DATADIR = '/Users/sebastiandresbach/data/s1Anfunco/Nifti'

# Set subjects to work on
subs = ['sub-07','sub-09','sub-10']
# Set sessions to work on
sessions = ['ses-01']


for sub in subs:

    # Define output folder
    outFolder = f'{DATADIR}/derivatives/{sub}/anat'
    # Create folder if not exists
    if not os.path.exists(outFolder):
        os.makedirs(outFolder)
        print("Subject directory is created")

    for ses in sessions:
        # Look for anatomical images.
        imgs = sorted(glob.glob(f'{DATADIR}/{sub}/{ses}/anat/{sub}_{ses}_highres-mp2rage_run-00*_uni.nii.gz'))

        # Load data of first run
        new = nb.load(imgs[0]).get_fdata()  # Load data
        affine = nb.load(imgs[0]).affine  # Load affine

        # We will register all images to the first
        # Read reference image in antsPy style
        fixed = ants.image_read(imgs[0])

        # Register images 2 and 3
        for n, img in enumerate(imgs[1:], start = 2):
            # Get a base name that we will use
            base = os.path.basename(img).rsplit('.', 2)[0]

            # Define moving image
            moving = ants.image_read(img)

            # Do registration
            mytx = ants.registration(fixed=fixed,
                                     moving=moving,
                                     type_of_transform='Rigid'
                                     )

            # Save image for quality control
            warped_moving = mytx['warpedmovout']
            ants.image_write(warped_moving, f'{outFolder}/{base}_registered.nii', ri=False)

            # Load registered images as array
            img = nb.load(f'{outFolder}/{base}_registered.nii').get_fdata()
            # Add it to new image
            new = np.add(new, img)

        # Divide new image by 3
        new = np.divide(new, 3)

        # Save averaged image
        ni_img = nb.Nifti1Image(new, affine)
        nb.save(ni_img,
                f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni.nii'
                )

        # Perform bias field correction
        img = ants.image_read(f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni.nii')
        img_n4 = ants.n4_bias_field_correction(img)

        # Save corrected image
        ants.image_write(img_n4,
                         f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected.nii',
                         ri=False
                         )

        # Because the N4 correction introduces some bright voxels, I will clip
        # the scale to 5000
        file = f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected.nii'
        nii = nb.load(f'{outFolder}/{sub}_{ses}_highres-mp2rage_average_uni_N4corrected.nii')
        affine = nii.affine
        header = nii.header
        data = nii.get_fdata()

        data[data > 5000] = 5000

        # Save clipped image
        ni_img = nb.Nifti1Image(data, affine=affine, header=header)
        nb.save(ni_img, f'{file}')
