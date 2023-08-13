function sigVec = gqcs_gaopin(dataX,snr,qcCoefs)
% Generate a quadratic chirp signal
% Pin Gao 2023.8.8

phaseVec = qcCoefs(1)*dataX + qcCoefs(2)*dataX.^2 + qcCoefs(3)*dataX.^3;
sigVec = sin(2*pi*phaseVec);
sigVec = snr*sigVec/norm(sigVec);
end

