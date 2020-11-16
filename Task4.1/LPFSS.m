clear;
% Script to simulate the step by step design of a lpf using a fss
% Step 1: Calculate the filter parameters: these are the value of the
%   lumped elements. Denormalize these filter parameters to the working
%   frequency.
% Step 2: 
%% Initial variables
path = '../../Task4/Report4_02/Images/';   % Path to save the files
resFile = 'results.txt';
file = fopen(resFile, 'w');

%% Low pass prototype
% Step 1: Set all the gk. The order of the filter is set to 5 using Pozar
% pag 407 because the order must be odd

N = 4;
g_lpnor = [1.6703 1.1926 2.3661 0.8419];

% 
% N = 5;
% g_lpnor = [1.7058 1.2296 2.5408 1.2296 1.7058];    % Pozar pag 406

printgk(g_lpnor, file, 'Normalized elements', 'g', 'F');

% All the elements are denormalized
w0 = 2*pi*7.5e9;
Z0 = 50;
g_lp = g_lpnor/(Z0*w0);

printgk(g_lp*1e12, file, 'Capacitors', 'C', 'pF');

%% Look for the patch lengths
load cvsl.mat
% The length is calculate using interpolation spline
Ldes = interp1(Copts, lengths, g_lp, 'spline');

% Prints and plots the length
plotcvsl(N, Copts*1e12, lengths, g_lp*1e12, Ldes, path, 'lengths');
printgk(Ldes, file, 'Patch Lengths', 'L', 'mm');
fclose(file);
