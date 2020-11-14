clear;
% Script to simulate the step by step design of a lpf using a fss
% Step 1: Calculate the filter parameters: these are the value of the
%   lumped elements. Denormalize these filter parameters to the working
%   frequency.
% Step 2: 
%% Initial variables
N = 3;
path = '../../Task4/Report4_02/Images/';   % Path to save the files
resFile = 'results.txt';
file = fopen(resFile, 'w');

%% Low pass prototype normalized
% Step 1: Set all the gk. The order of the filter is set to 4 using Pozar
% pag 407
N = 4;
g_lpnor = [1.6703 1.1926 2.3661 0.8419 1.9841];    % Pozar pag 406

% Step 2:
w = 0:0.01:(2*pi);  % Normalized frequency
F = calculateFLPF(N, w, g_lpnor);

% Step 3
[s11, s21] = multFretS(N, length(w), F);

% Step 4
plotNormalized(w, s11, s21, path, 'LowPass', 'Normalized f [Hz]');
printgk(g_lpnor, file, 'Low Pass prototype Normalized');

%% Denormalize low pass filter
fc = 7.5e9; % Cutoff frequency (Hz)
w0 = 2*pi*fc;
w = 2*pi*(0:10e6:15e9);

Z0 = 1;

g_lp = zeros(1, length(g_lpnor));
g_lp(1:2:end-1) = g_lpnor(1:2:end-1)/(Z0*w0);    % Denormalizing capacitors
g_lp(2:2:end-1) = g_lpnor(2:2:end-1)*Z0/w0;  % Denormalizing inductances
g_lp(end) = g_lpnor(end)/Z0;                     % Denormalizing the final resistance

F = calculateFLPF(N, w, g_lp);
[s11, s21] = multFretS(N, length(w), F);

plotNormalized(w/1e9, s11, s21, path, 'LowPass', 'f [GHz]');
printgk(g_lp, file, 'Low Pass prototype Normalized');

fclose(file);