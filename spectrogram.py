
import matplotlib.pyplot as plt
import numpy as np
from read_plot import *

t, w, amp, v = read_plot("omega.data")

s1 = []
for i in range(len(t)) :
    s1.append(amp[i]*np.sin(w[i]*t[i]))
    v[i] /= 3.e8

x = s1
Fs = 4096
NFFT = int(Fs/16.0)
NOVL = int(NFFT*15/16.0)

fig, (ax1, ax2) = plt.subplots(nrows=2)
ax1.plot(t,x, linewidth=.7)
ax2.plot(t,v, linewidth=.7)
#Pxx, freqs, bins, im = ax2.specgram(x, NFFT=NFFT, Fs=Fs, noverlap=NOVL)
ax1.set_ylabel('Strain')
ax1.set_xlabel('Time (s)')
ax2.set_ylabel('v/c')
ax2.set_xlabel('Time (s)')
#fig.colorbar(im)
plt.show()


Pxx, freqs, bins, im = plt.specgram(x,NFFT=NFFT, Fs=Fs, noverlap=NOVL, xextent=[0,5])
plt.ylabel('Frequency')
plt.xlabel('Samples (~2*t)')
plt.colorbar(im)
plt.show()
