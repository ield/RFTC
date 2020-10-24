%% Author: ield
% Script created to calculate the dispersion parameters of a given
% waveguide of length a, b (a>b), in the temn and tmmn modes, between 
% fmin and fmax
function [alpha,beta, fc, kcut] = dispersion(m, n, a, b, c, fmin, fmax)

f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)

kc = sqrt((m*pi/a)^2+(n*pi/b)^2);   % Formula pag 17 pozar
fc = kc*c/(2*pi) ;               % Formula 3.84 Pozar
kcut = 2*pi*fc/c;

k = 2*pi*f/c;                       % Formula pag 17 pozar
beta = sqrt(k.^2 - kc^2);           % Formula pag 17 pozar
alpha = sqrt(kc^2-(2*pi*f).^2./c^2);% Formula slide 26 moodle presentation 3.1
end

