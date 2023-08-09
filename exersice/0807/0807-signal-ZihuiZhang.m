%Signal parameters
qcCoefs = [10,3,3];
A = 10;

samplIntrvl = 1e-2;
timeVec = 0:samplIntrvl:1.0;

phaseVec = qcCoefs(1)*timeVec + qcCoefs(2)*timeVec.^2 + qcCoefs(3)*timeVec.^3;
sigVec = A*sin(2*pi*phaseVec);

figure
plot(timeVec,sigVec,'Marker','.','MarkerSize',10);
xlabel('Time (sec)');
title('Sampled signal');