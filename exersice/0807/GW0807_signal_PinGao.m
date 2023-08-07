% generate a analogous GW signal
% PinGao 2023.8.7

clear
a1=10;
a2=3;
a3=3;
A = 10;
samplIntrvl = 0.001; %seconds
sigLen = 1.0; %seconds
timeVec = 0:samplIntrvl:sigLen;
phaseVec = a1*timeVec + a2*timeVec.^2 + a3*timeVec.^3;
sigVec = A * sin(2*pi*phaseVec);
figure
plot(timeVec,sigVec,'b.');
xlabel('T')
ylabel('A')
