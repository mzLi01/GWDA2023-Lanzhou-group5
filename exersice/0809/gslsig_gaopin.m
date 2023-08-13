% Use the sinusoidal signal to generate simulated LIGO signal
% gaopin 2023.8.9

clear
%% F_+ and F_x calculation
%Polar
theta = 0.5; %Select signal Angle (radian)
%Azimuthal
phi = 0.5;

[A,D] = meshgrid(phi,theta); %极坐标系转换为直角坐标系
X = sin(D).*cos(A);
Y = sin(D).*sin(A);
Z = cos(D);

%Generate function values
fPlus = zeros(length(theta),length(phi));
fCross = zeros(length(theta),length(phi));
for lp1 = 1:length(phi)
    for lp2 = 1:length(theta)
        [fPlus(lp2,lp1),fCross(lp2,lp1)] = detframefpfc(theta(lp2),phi(lp1));
    end
end
%% generate LIGO signal
% signal in Wave frame conventions
a = 1;
b = 1;
f = 20;
t = 0:0.01:1;
phi = 0.5;
h1 = a*sin(2*pi*f*t);
h2 = b*sin(2*pi*f*t+phi);

% signal in detector frame conventions
h11 = h1*fPlus;
h22 = h2*fCross;

%Plot
figure;
plot(h11)
figure 
plot(h22)
% surf(X,Y,Z,abs(fPlus));
% shading interp;
% axis equal;
% colorbar;
% figure;
% surf(X,Y,Z,abs(fCross));
% shading interp;
% axis equal;
% colorbar;