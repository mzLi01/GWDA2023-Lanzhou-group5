import numpy as np

def crcbgenqcsig(dataX,snr,qcCoefs):
# Generate a quadratic chirp signal
# S = CRCBGENQSIG(X,SNR,C)
# Generates a quadratic chirp signal S. X is the vector of
# time stamps at which the samples of the signal are to be computed. SNR is
# the matched filtering signal-to-noise ratio of S and C is the vector of
# three coefficients [a1, a2, a3] that parametrize the phase of the signal:
# a1*t+a2*t^2+a3*t^3. 

#Soumya D. Mohanty, May 2018

    phaseVec = qcCoefs[0]*dataX + qcCoefs[1]*dataX**2 + qcCoefs[2]*dataX**3
    sigVec = np.sin(2*np.pi*phaseVec)
    #print(sigVec)
   # print(np.sqrt(sigVec**2))
    sigVec = snr*sigVec/np.linalg.norm(sigVec, 2)
    return sigVec