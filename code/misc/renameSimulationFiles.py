import subprocess
import glob

ROOT = '/Users/sebastiandresbach/data/s1Anfunco/Nifti/derivatives/designFiles/stimulationPatterns'

subs = ['sub-15', 'sub-16', 'sub-17', 'sub-18']

for sub in subs:
    subFolder = f'{ROOT}/{sub}'

    files = sorted(glob.glob(f'{subFolder}/*ses-00*_D*.txt'))

    for file in files:
        newName = file[:-28] + file[-27:]
        command = f'cp {file} {newName}'
        subprocess.run(command, shell=True)
