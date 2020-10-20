clear;
%% Author: ield
% Script created to solve the Exercise 1 of task 2 of the RFTC subject of
% the MSTC

%% a
% Constant declaration
a = 9.05e-3;        % a [m]
b = 9.525e-3;       % b [m]
e = 1;              % Dielectric constant
c = 3e8;            % Speed of light [m/s]
fmin = 1000;
fmax = 3000;
f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)

% Plots
figure('Color',[1 1 1]);
ylim([f(1)/1e9 f(end)/1e9]);
xlabel('\alpha, \beta');
ylabel('Frequency [GHz]');
hold on;

% First mode: according to Pozar TE10
[alpha, beta] = dispersion(1, 0, a, b, c, fmin, fmax);
plot(beta, f/1e9);
plot(alpha, f/1e9);
xlim([0 max(beta)]);

% Second mode: according to Pozar: TM11
[alpha, beta] = dispersion(1, 1, a, b, c, fmin, fmax);
plot(beta, f/1e9);
plot(alpha, f/1e9);

%%