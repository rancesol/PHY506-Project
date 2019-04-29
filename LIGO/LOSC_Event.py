
eventname = ''
eventname = 'GW150914'
plottype = 'pdf'

# Standard python imports:
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
from scipy.signal import butter, filtfilt, iirdesign, zpk2tf, freqz
import h5py
import json
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# LIGO-specific import:
import readligo as rl

# Read event properties from a local json file
fnjson = 'BBH_events_v3.json'
try:
    events = json.load(open(fnjson,'r'))
except IOError:
    print('Cannot find resource file ' + fnjson)
    print('Quitting.')
    quit()

# Check for choice of eventname:
try:
    events[eventname]
except:
    print('You must select an eventname that is in ' + fnjson + '! Quitting.')
    quit()

# Extract parameters from event:
event = events[eventname]
fn_H1 = event['fn_H1']
fn_L1 = event['fn_L1']
fn_template = event['fn_template']
fs = event['fs']
tevent = event['tevent']
fband = event['fband']
print('Reading in parameters for event ' + event['name'])
#print(event)

# Read in event data
try:
    strain_H1, time_H1, chan_dict_H1 = rl.loaddata(fn_H1, 'H1')
    strain_L1, time_L1, chan_dict_L1 = rl.loaddata(fn_L1, 'L1')
except:
    print('cannot find data files!')
    print('Quitting.')
    quit()

# H1 and L1 have the same time vector
time = time_H1
# time interval with uniform sampling
dt = time[1] - time[0]

# select the a time window around the event
deltat = 5
indxt = np.where((time >= tevent-deltat) & (time < tevent+deltat))

make_psds = 1
if make_psds:
    # number of sample for the FFT:
    NFFT = 4*fs
    Pxx_H1, freqs = mlab.psd(strain_H1, Fs = fs, NFFT=NFFT)
    Pxx_L1, freqs = mlab.psd(strain_L1, Fs = fs, NFFT=NFFT)

    # interpolations of ASDs computed for whitening
    psd_H1 = interp1d(freqs, Pxx_H1)
    psd_L1 = interp1d(freqs, Pxx_L1)

# The signal is dominated by low frequency noise.
# We need to do some signal processing.


# This function will whiten the data
def whiten(strain, interp_psd, dt):
    Nt = len(strain)
    freqs = np.fft.rfftfreq(Nt, dt)
    freqs1 = np.linspace(0,2048.,Nt/2+1)

    hf = np.fft.rfft(strain)
    norm = 1./np.sqrt(1./(dt*2))
    white_hf = hf / np.sqrt(interp_psd(freqs)) * norm
    white_ht = np.fft.irfft(white_hf, n=Nt)
    return white_ht

whiten_data = 1
if whiten_data:
    strain_H1_whiten = whiten(strain_H1,psd_H1,dt)
    strain_L1_whiten = whiten(strain_L1,psd_L1,dt)

    bb, ab = butter(4, [fband[0]*2./fs, fband[1]*2./fs], btype='band')
    normalization = np.sqrt((fband[1]-fband[0])/(fs/2))
    strain_H1_whitenbp = filtfilt(bb, ab, strain_H1_whiten) / normalization
    strain_L1_whitenbp = filtfilt(bb, ab, strain_L1_whiten) / normalization


# pick shorter FFT time interval
NFFT = int(fs/16.0)
# with a lot of overlap for short-time features
NOVL = int(NFFT*15./16.)
# choose window that minimizes 'spectral leakage'
window = np.blackman(NFFT)

# color map
spec_cmap = 'ocean'

# plot H1 spectrogram
plt.figure(figsize=(10,6))
spec_H1, freqs, bins, im = plt.specgram(strain_H1_whiten[indxt], NFFT=NFFT, Fs=fs, window=window,
        noverlap=NOVL, cmap=spec_cmap, xextent=[-deltat,deltat])
plt.xlabel('time (s) since ' + str(tevent))
plt.ylabel('Frequency (Hz)')
plt.colorbar()
plt.axis([-0.5, 0.5, 0, 500])
plt.title('aLIGO H1 strain data near ' + eventname)
plt.show()

