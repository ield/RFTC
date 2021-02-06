% Engineer: ield
% Script to complete task 9
% Constants in the radar equation
clear;
sigma = 1;                      % Radar cross section
T = 9;                          % Scanning time
KT0 = 4e-21;                    % Boltzmann and reference temperature
F = 2.5;  nf = 10*log10(F);        % Noise figure
L = 8;  l = 10*log10(L);        % Losses
SNR = 13;   snr = 10*log10(SNR);     
Omega = 2*pi/9;                 % Solid angle (sr) reasoned in notes RFTC-12
eff = 0.5;                      % Efficiency of the radiating element
dc = 0.2;                       % Duty cycle
c = 3e8;                        % Speed of light (m/s)
R = 100e3;                      % Range (m)

path = '../../Task9_PhasedArray/Images/';
%% Number of elements
theta_a = 8;        %º
theta_a = theta_a*pi/180;

Nel_aperture = ceil(4/(theta_a^2 * eff));
fprintf('\nThere are needed at least %i elements per section\n', ...
        Nel_aperture);

% The aperture must be square so
Nel_size = ceil(sqrt(Nel_aperture));
Nel_aperture = Nel_size^2;
fprintf('They are grouped %i elements in sqaures of %i x %i\n', ...
        Nel_aperture, Nel_size, Nel_size);

%% Calculate the working frequency
f = linspace(500e6, 10e9, 5000);     % Frequency aquis (Hz)
lambda = c./f;                      % Wavelength (m)

CTE = sigma*T./(4*pi*KT0*nf*l*snr*Omega)*eff/4*dc;

Pel = R^4./(lambda.^2*Nel_aperture.*CTE);    % Power needed at each element W
Pel = 10*log10(Pel*1000);                    % Power needed at each element (dBm)

area = Nel_aperture*lambda.^2/4;

set(0, 'DefaultAxesFontName', 'Times New Roman');
figure('Color',[1 1 1]);
plot(f/1e6, Pel, 'r');

ylabel('Power (dBm)');
xlabel('Frequency (MHz)');
xlim([f(1) f(end)]/1e6);

yyaxis right
plot(f/1e6, area, 'b'); hold on;
plot(f/1e6, ones(size(f)), 'b--');
ylabel('Area (m^2)');

legend('Power', 'Area', 'Location', 'north');
saveas(gca, [path, 'power_area'],'epsc');

pos_work = find(area < 1, 1);
work_f = f(pos_work);
fprintf('The working frequency is %f MHz\n', work_f/1e6);

%% Fixing the area, selecting the minimum power
area = 1;

theta_a = linspace(0.00001, 8, 5000);        %º
theta_a = theta_a*pi/180;

Nel_aperture = ceil(4./(theta_a.^2 * eff));
Nel_size = ceil(sqrt(Nel_aperture));
Nel_aperture = Nel_size.^2;

lambda = sqrt(4*area./Nel_aperture);
f = c./lambda;

CTE = sigma*T./(4*pi*KT0*nf*l*snr*Omega)*eff/4*dc;
Pel = R^4./(lambda.^2.*Nel_aperture.*CTE);
Pel = 10*log10(Pel*1000); 

figure('Color',[1 1 1]);
plot(theta_a*180/pi, Pel, 'k');

ylabel('Power (dBm)');
xlabel('\theta_a = \theta_h (º)');
xlim([theta_a(1) theta_a(end)]*180/pi);
ylim([0 100]);

% saveas(gca, [path, 'power_theta'],'epsc');

fprintf('Each element has a power of %f dBm\n', mean(Pel));
%% Phase shifters
theta_a = 8;        %º
theta_a = theta_a*pi/180;

Nel_aperture = ceil(4/(theta_a^2 * eff));
fprintf('\nThere are needed at least %i elements per section\n', ...
        Nel_aperture);

% The aperture must be square so
Nel_size = ceil(sqrt(Nel_aperture));
Nel_aperture = Nel_size^2;

% Number of explorations
explorations_needed = ceil(Omega*Nel_aperture*eff*pi/8);
sr_covered_exploration = Omega / explorations_needed;
fprintf('There are done %i explorations convering %f sr in each exploration\n', ...
        explorations_needed, sr_covered_exploration);
    
% Phase shifter bits
N_bits = ceil(log2(explorations_needed));
fprintf('For %i explorations, they are needed phase shifters of %i bits\n', ...
        explorations_needed, N_bits);

%% T/R module
% Run after Calculate the working frequency
p_max = 53;         % Maximum power (dBm)

pos_work = find(Pel > p_max, 1);
work_f = f(pos_work);
fprintf('f = %f MHz; A = %f m^2\n', work_f/1e6, area(pos_work));


%%
F_target = 2;           f_target = 10^(F_target/10);

F_circ = 1;               f_circ = 10^(F_circ/10);

F_lim = 0.7;              f_lim = 10^(F_lim/10);

f_max_lna = (f_target - f_circ - (f_lim - 1)*f_circ)/(f_circ*f_lim) + 1;
F_max_lna = 10*log10(f_max_lna);
fprintf('F max LNA = %f dB\n', F_max_lna);

%% Complete noise figure
F_lna = 0.28;
G_lna = 19;

F_bpf = 3.2;

F_bga2866 = 1.1;
G_bga2866 = 23.9;

F_swi = 0.35;

F_ps = 4;

% They are put all the F and G together: the gain of the attenuators is its
% F with - sign. 
F_all = [F_circ F_lim F_lna F_bpf F_bga2866 F_bpf F_swi F_ps];
G_all = [-F_circ -F_lim G_lna -F_bpf G_bga2866 -F_bpf F_swi -F_ps];

% All the f and g are transformed to linear units
F_all = 10.^(F_all/10);
G_all = 10.^(G_all/10);

% The first element in G is set to 1 so the multiplication is done correctly
% G_all = [1 G_all];

F = F_all(1);          % Global noise figure

for ii = 2:length(F_all)
    F = F + (F_all(ii) - 1)/prod(G_all(1:ii-1));    
end

F = 10*log10(F);

fprintf('F = %f dB\n', F);










