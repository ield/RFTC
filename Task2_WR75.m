clear;
%% Author: ield
% Script created to solve the Exercise 1 of task 2 of the RFTC subject of
% the MSTC

%% a
% Constant declaration
a = 19.05e-3;        % a [m]
b = 9.525e-3;       % b [m]
e = 1;              % Dielectric constant
c = 3e8;            % Speed of light [m/s]
fmin = 0;
fmax = 2000;
f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)

% Plots
figure('Color',[1 1 1]);
ylim([f(1)/1e9 f(end)/1e9]);
xlabel('\alpha, \beta');
ylabel('Frequency [GHz]');

hold on;

% First mode: according to Pozar TE10
[alpha1, beta1, fc1, kcut1] = dispersion(1, 0, a, b, c, fmin, fmax);
plot(beta1, f/1e9);
plot(alpha1, f/1e9);
xlim([0 max(beta1)]);

% Second mode: according to Pozar: TM11
[alpha2, beta2, fc2, kcut2] = dispersion(2, 0, a, b, c, fmin, fmax);
plot(beta2, f/1e9);
plot(alpha2, f/1e9);

legend('\beta TE_1_0', '\alpha TE_1_0', '\beta TE_0_1', '\alpha TE_0_1');
saveas(gcf, '../Task2/Images/Ex1a', 'svg');
hold off;

%% Ex1.b
% Remake the graph and plot the special points
figure('Color',[1 1 1]);
xlim([0 max(beta1)]);
ylim([f(1)/1e9 f(end)/1e9]);
xlabel('\alpha, \beta');
ylabel('Frequency [GHz]');
hold on;

plot(beta1, f/1e9, ':');
plot(alpha1, f/1e9, ':');
plot(beta2, f/1e9, ':');
plot(alpha2, f/1e9, ':');
legend('\alpha TE_1_0', '\beta TE_1_0', '\alpha TE_0_1', '\beta TE_0_1')
hold on;

plot(0, fc1/1e9, 'o', 'LineWidth', 2)
plot(0, fc2/1e9, 'o', 'LineWidth', 2)
plot(kcut1, 0, 'o', 'LineWidth', 2)
plot(kcut2, 0, 'o', 'LineWidth', 2)

hold on;
smbw = [fc1 fc2]/1e9
plot([0 0], smbw, 'LineWidth', 2);

recbw = [1.25*fc1 0.95*fc2]/1e9;
plot([0 0], recbw, 'LineWidth', 4);

plot(0, 1.25*fc1/1e9, 'o', 'LineWidth', 2)
plot(0, 0.95*fc2/1e9, 'o', 'LineWidth', 2)

legend('\beta TE_1_0', '\alpha TE_1_0', '\beta TE_0_1', '\alpha TE_0_1',...
    'f_cTE_1_0', 'f_cTE_0_1', 'k_{cut}TE_1_0', 'k_{cut}TE_0_1',...
    'Single-Mode Bandwidth', 'Recommended Bandwidth', 'f_{min}', 'f_{max}');
saveas(gcf, '../Task2/Images/Ex1b', 'svg');

%% Ex1.c
fmax = 0.95*fc2;
k_fmax = 2*pi*fmax/c;
kc = sqrt((1*pi/a)^2+(0*pi/b)^2);
betamax = sqrt(k_fmax^2 - kc^2);
vp = 2*pi*fmax / betamax

Vg_1 = c*sqrt(1-(fc1/fmax)^2)

fmax_ = fmax - 10e6;
k_fmax_ = 2*pi*fmax/c;
kc = sqrt((1*pi/a)^2+(0*pi/b)^2);
betamax = sqrt(k_fmax^2 - kc^2);

%% Ex1.d
% The transmitted power is obtained from diapo 56 in slides 2.1
% Constant declaration
a = 19.05e-3;        % a [m]
b = 9.525e-3;       % b [m]
e = 1;              % Dielectric constant
c = 3e8;            % Speed of light [m/s]
eta = 120*pi;
fmin = 0;
fmax = 2000;
f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)

% Plots
figure('Color',[1 1 1]);
xlim([f(1)/1e9 f(end)/1e9]);
xlabel('Frequency [GHz]');
ylabel('Transmitted Power [MW]');

hold on;

% Obtaining parameters for equation:
% gamma
[alpha1, beta1, fc1, kcut1] = dispersion(1, 0, a, b, c, fmin, fmax);
m = 1;
n = 0;
gamma = real(alpha1) + 1i*real(beta1);

% kc
kc = sqrt((m*pi/a)^2+(n*pi/b)^2);

% ZTe
Zte = obtainZte(kc, f);

% Obtain integral
Emax = 3e6;
integral = intHHcong10(a, b, Emax);

%Solve equiation
W = (gamma .* conj(gamma))/(2*kc^2).*real(Zte)*integral;

plot(f/1e9, W/1e6, 'k', 'LineWidth', 1);


% Find the transmited power at fmax
fmax_u = 0.95*fc2;
index = find(f' > fmax_u);
plot(f(index(1))/1e9, W(index(1))/1e6, 'o', 'LineWidth', 2);
legend('Transmitted power TE_1_0', 'Transmitted power in f_{max}');
saveas(gcf, '../Task2/Images/Ex1e', 'svg');

%% Ex 1.f
a = 19.05e-3;        % a [m]
b = 9.525e-3;       % b [m]
er = 1;              % Dielectric constant
c = 3e8;            % Speed of light [m/s]
eta = 120*pi;
fmin = 0;
fmax = 2000;
f = 10e6*[fmin:fmax];     % Frequency axis [Hz] (from 0 to 10 GHz with jump of 10 MHz)

% Obtain alphad
cond = 5.8e7;
e0 = 8.854e-12;           % Constant. See pozar appendix

% % epr = e0*er;              % Pozar pag10
% epr = er;
% tand = cond./(2*pi*f*epr);% Pozar pag10
% 
% k = 2*pi*f/c;
% [~, beta1, fc, ~] = dispersion(1, 0, a, b, c, fmin, fmax);
% alphad = k.^2.*tand./beta1;     % Pozar pag 117
% 
% % Obtain alphac
mu0 = 4*pi*1e-7;
Rs = sqrt((2*pi*f*mu0)/(2*cond));% Pozar pag28. Rs = 1/(sigmadelta)

alphac = 1./eta.*Rs.*((fc./f).^2+a/(2*b))./((a/2).*sqrt(1-(fc./f).^2)); %Microwave T3 slides.


alphaTot = alphac *8.685889638; % Change to dB

figure('Color',[1 1 1]);

plot(f/1e9, alphaTot, 'k', 'LineWidth', 1);
xlim([f(1)/1e9 f(end)/1e9]);
xlabel('Frequency [GHz]');
ylabel('Attenuation constant [dB/m]');
saveas(gcf, '../Task2/Images/Ex1f', 'svg');

