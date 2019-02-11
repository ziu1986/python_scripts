import numpy as np
from scipy import fftpack
import matplotlib.pyplot as plt

plt.close('all')

f = 10  # Frequency, in cycles per second, or Hertz
f_s = 200  # Sampling rate, or number of measurements per second

t = np.linspace(0, 2, 2 * f_s, endpoint=False)
x = np.sin(f * 2 * np.pi * t)+np.cos((f*0.3) * 2 * np.pi * t)+np.sin((f*4.4) * 2 * np.pi * t)
x1 = np.sin(f * 2 * np.pi * t)

fig2 = plt.figure(2)
ax21 = plt.subplot(211)
ax21.plot(t, x)
ax22 = plt.subplot(212)
ax22.plot(t, x1)
ax22.set_xlabel('Time [s]')
ax22.set_ylabel('Signal amplitude', y=1);

X = fftpack.fft(x)
freqs = fftpack.fftfreq(len(x)) * f_s

X1 = fftpack.fft(x1)
freqs1 = fftpack.fftfreq(len(x1)) * f_s

fig1 = plt.figure(1)
ax11 = plt.subplot(211)
ax11.stem(freqs, np.abs(X))
ax12 = plt.subplot(212)
ax12.stem(freqs1, np.abs(X1))
ax12.set_xlabel('Frequency in Hertz [Hz]')
ax12.set_ylabel('Frequency Domain (Spectrum) Magnitude', y=1)
for ax in fig1.axes:
    ax.set_xlim(0, f_s / 2)
    ax.set_ylim(-5, 110)

plt.show(block=False)
