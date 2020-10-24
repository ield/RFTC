%% Author: ield
% Script created to solve the Exercise 2 of task 2 of the RFTC subject of
% the MSTC

%% Calculate the transmitted power
a = 20.24e-3;        % a [m]
er = 1;              % Dielectric constant
c = 3e8;            % Speed of light [m/s]
eta = 120*pi;
fmin = 0;
fmax = 1500;
f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)
E0 = 1.5e6;
p11 = 3.832

% Calculate transmitted power
fc = c/(2*pi)*p11/a
W = 0.995*pi*a^2/(4*eta)*E0^2.*sqrt(1-(fc./f).^2);


figure('Color',[1 1 1]);

plot(f/1e9, W/1e6, 'k', 'LineWidth', 1);

xlim([f(1)/1e9 f(end)/1e9]);
xlabel('Frequency [GHz]');
ylabel('Transmitted Power [MW]');

saveas(gcf, '../Task2/Images/Ex2a.eps');

%% Calculate attenuation constant
fmax = 2000;
f = 10e6*[fmin:fmax];

% Obtain alphad
cond = 6.17e7;
e0 = 8.854e-12;           % Constant. See pozar appendix
mu0 = 4*pi*1e-7;


Rs = sqrt((2*pi*f*mu0)/(2*cond));% Pozar pag28. Rs = 1/(sigmadelta)
alphac = Rs./eta.*((fc./f).^2+1/(p11^2-1))./(a*sqrt(1-(fc./f).^2)); %Microwave T3 slides.


alphaTot = alphac *8.685889638; % Change to dB

figure('Color',[1 1 1]);

plot(f/1e9, alphaTot, 'k', 'LineWidth', 1);
xlim([f(1)/1e9 f(end)/1e9]);
xlabel('Frequency [GHz]');
ylabel('Attenuation constant [dB/m]');
saveas(gcf, '../Task2/Images/Ex2b.eps');
