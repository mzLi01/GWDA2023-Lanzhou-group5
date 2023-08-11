import cryptography
import numpy as np
import matplotlib.pyplot as plt
#Plot the quadratic chirp signal
#Signal parameters
a1=10
a2=3
a3=3
A = 10
#Instantaneous frequency after 1 sec is 
maxFreq = a1+2*a2+3*a3
#Nyqust frequency guess: 2 * max. instantaneous frequency
nyqFreq = 2*maxFreq
#Sampling frequency
samplFreq = 5*nyqFreq
samplIntrvl = 1/samplFreq

# Time samples
timeVec = np.arange(0,1/samplIntrvl)*samplIntrvl
# Number of samples
nSamples = len(timeVec)

# Generate the signal
sigVec = crcbgenqcsig(timeVec,A,[a1,a2,a3]);

plt.plot(timeVec,sigVec)
plt.xlabel('Time (sec)')

plt.title('Sampled signal')
#plt.savefig("crcbgenqcsig_example.pdf")