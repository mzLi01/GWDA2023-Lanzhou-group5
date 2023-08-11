#Signals and audio
# Time domain signal: sum of two sinusoids with frequencies 5 Hz and 1 Hz
nSamples = 20480
samplingFreq = 100
#t = (0:(nSamples-1))/samplingFreq
t = np.arange(nSamples-1)/samplingFreq
def my_signal_2(t):
    
    st = np.sin(2*np.pi*2*t)*(np.sin(2*np.pi*5*t)+2*np.sin(2*5.5*np.pi*t))
    return 