clear all
close all
clc

epsinf = 3.77;
gamma = 0.05;
wp = 9.1;

omega = linspace(1,6,251);

epsilon = epsinf - (wp^2)./(omega.^2 + 1i.*gamma.*omega);
epsr = real(epsilon);
epsi = imag(epsilon);
epsmag = sqrt(epsr.^2 + epsi.^2);

n = sqrt((epsmag + epsr)/2);
k = sqrt((epsmag - epsr)/2);

plot(omega,epsr,omega,epsi)

thing_to_save = [omega',n',k'];
save('silver_normal.dat','thing_to_save','-ascii')
