""" Plotting stimulation sound for manuscript """

from scipy.io import wavfile
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('dark_background')

soundFile = "30_5_unpredictable.wav"
samplerate, data = wavfile.read(soundFile)

x = np.linspace(0, len(data), len(data))

fig = plt.figure(figsize=(100, 2))
plt.plot(x, data, color='white')
plt.axis('off')
plt.savefig(f'plots/stimulationSound.png', bbox_inches='tight', pad_inches=0)
plt.show()

# Plot zoomed
fig = plt.figure(figsize=(100, 2))
plt.plot(x[:int(len(data)/6)], data[:int(len(data)/6)], color='white', linewidth=3)
plt.axis('off')
plt.savefig(f'plots/stimulationSoundZoomed.png', bbox_inches='tight', pad_inches=0)
plt.show()