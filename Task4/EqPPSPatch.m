%% Engineer: ield

clear;
path = '../../Task4/Images';    % Path to save images
%% Step 1. Load in allData all the measures taken.
% all data is an array of structures that each element
%   l: length of the patch
%   f: frequencies evaluated
%   s11: s11 measures (real and imag)
%   s31: s31 measures (real and imag)
allData = textToSignal('s11_all.txt', 's31_all.txt');
Z0 = 50;

%% Step 2. Transform the s parameters to an admittance
% The F par of an admittance are found in pozar pag 190
% The value of the admittance is C, and C in terms of s par is found in
% pozar pag 192
% Since the quadripole is symmetric for v and h, the s31 is the s21 if we
% were considering a quadripole.
% Since the quadripole is symmetric, s11 = s22 and s12 = s21
Y = 1/Z0 * ((1-allData(1).s11).^2-allData(1).s31.^2)./(2*allData(1).s31);
plotAdmittance(allData(1).f, Y, 'Y(f)', path, 'y(f)')





