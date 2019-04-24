
import matplotlib.pyplot as plt
import numpy as np
from read_plot import *

np.random.seed(19680801)

t, w, amp = read_plot("omega.data")

dt = 0.01
s1 = []
for i in range(len(t)) :
    s1.append(amp[i]*np.sin(w[i]*t[i]))

nse = 0.01*np.random.random(size=len(t))

x = s1 + nse
NFFT = 1024
Fs = int(5.0 / dt)

fig, (ax1, ax2) = plt.subplots(nrows=2)
ax1.plot(t,x, linewidth=.7)
Pxx, freqs, bins, im = ax2.specgram(x, NFFT=NFFT, Fs=Fs, noverlap=1000)
plt.show()
