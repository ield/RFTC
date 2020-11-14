%% Engineer: ield

clear;
path = '../../Task4/Images/';    % Path to save images
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

%% Step 3. Compare to a parallel capacitor in 0-7 GHz
f7g = allData(1).f(1:find(allData(1).f >= 7, 1))*1e9;   % Frequencies from 0-7GHz
Y7g = Y(1:find(allData(1).f >= 7, 1));              % Y values from 0-7GHz scaled to use fminunc

%It optimzes the value of C so that it adjusts the desired function
options=optimoptions('lsqnonlin');
% options.Display='iter';

Cini = 1;                          % Initial value of C
clc
[Copt] = lsqnonlin(@(C) minDisC(C, f7g, imag(Y7g)), Cini, [], [], options);
Copt = Copt*1e-12;

Y_copt = 1i*2*pi*allData(1).f*1e9*Copt;
plotAdmittanceAndC(allData(1).f, Y, Y_copt, '', path, 'y_c(f)');
fprintf('C = %f pF', Copt*1e12);

%% Step 3. Compare to a parallel capacitor + inductance
f = allData(1).f*1e9;   % All frequencies
Y_des = Y;              % Y values scaled to use fminunc

%It optimzes the value of C and L so that it adjusts the desired function
LCini = [1 1];
clc
LCopt = lsqnonlin(@(LC) minDisLC(LC, f, imag(Y_des)), LCini, [], [], options);
L = LCopt(1)*1e-9;
C = LCopt(2)*1e-12;
Y_copt = 1i*2*pi*f*C./(1-(2*pi*f).^2*L*C);
plotAdmittanceAndC(allData(1).f, Y, Y_copt, '', path, 'y_lc(f)');
fprintf('L = %f nH and C = %f pF', L*1e9, C*1e12);

%% Step 4. Model of a single C for all lengths
% Represent all the values of optimum C
f7g = allData(1).f(1:find(allData(1).f >= 7, 1))*1e9;   % Frequencies from 0-7GHz
%It optimzes the value of C so that it adjusts the desired function
Cini = 1;                          % Initial value of C

lengths = zeros(1, length(allData)-1);
Copts = zeros(1, length(allData)-1);
for ii = 1:length(allData)-1
    lengths(ii) = allData(ii).length;
    Y = 1/Z0 * ((1-allData(ii).s11).^2+allData(ii).s31.^2)./(2*allData(ii).s31);
    Y7g = Y(1:find(allData(ii).f >= 7, 1));             % Y values from 0-7GHz scaled to use fminunc

    
    [Copt] = lsqnonlin(@(C) minDisC(C, f7g, imag(Y7g)), Cini, [], [], options);
    Copt*1e-12;
    Copts(ii) = Copt*1e-12;
end

plotcvsl(lengths, Copts*1e12, 'Optimum C (pF)', path, 'cvsl');









