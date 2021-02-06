% Engineer: ield
% Script to complete task 9
%% Question 1
% Constants in the radar equation
sigma = 1;                      % Radar cross section
T = 9;                          % Scanning time
KT0 = 4e-21;                    % Boltzmann and reference temperature
F = 3;  f = 10*log10(F);        % Noise figure
L = 8;  l = 10*log10(L);        % Losses
SNR = 13;   snr = 10*log10(SNR);     
Omega = pi^2/9;                 % Solid angle (sr) reasoned in notes RFTC-12
eff = 0.5;                      % Efficiency of the radiating element
dc = 0.2;                       % Duty cycle

CTE = sigma*T/(4*pi*KT0*f*l*snr*Omega)*eff/4*dc;

brands = {'TGA2517', 'TIM3742', 'FLL410IK', 'TGA2514'};

f = [10 4 2.6 15]*1e9;          % Frequency for each brand (Hz)
c = 3e8;                        % Speed of light
lambda = c./f;                  % Wavelength for each brand (m)

Pel = [42 45 46 38];            % Power of each brand (dBm)
pel = 10.^((Pel-30)/10);        % Power of each brand (W)

R = [100 200 300 400 500]*1e3;  % Range (m)

for ii = 1:length(R)
    Nel = 3*ceil(R(ii).^2./(lambda.*sqrt(pel*CTE)));    % Round up because you need integer elements
    
    fprintf('\n\\item For a range of %i km there are needed: ', R(ii)/1e3);
    for jj =1:length(brands)
        fprintf('%i %s, ', Nel(jj), brands{jj});
    end
end

%% Question 2
brand_number = 3;           % It is selected TIM3742, for 500 km
Nel_aperture = ceil(R(end).^2./(lambda(brand_number).*sqrt(pel(brand_number)*CTE)));
fprintf('\nThere are needed at least %i elements of brand %s per section\n', ...
        Nel_aperture, brands{brand_number});

% The aperture must be square so
Nel_size = ceil(sqrt(Nel_aperture));
Nel_aperture = Nel_size^2;
fprintf('There are grouped %i elements in sqaures of %i x %i\n', ...
        Nel_aperture, Nel_size, Nel_size);
    
% Beamwidth
bw = sqrt(4/(Nel_aperture*eff));
fprintf('bw_e = %f rad. bw_h = %f rad\n', bw, bw);

% Number of explorations
elevation_angle = pi/6;       % Elevation angle (rad)
azimuth_angle = 2*pi/3;       % Azimuth angle (rad)
explorations_elevation = ceil(elevation_angle/bw);
explorations_azimuth = ceil(azimuth_angle/bw);
explorations_needed = explorations_elevation * explorations_azimuth;
% explorations_needed = ceil(solid_angle*Nel_aperture*eff*pi/8);
sr_covered_exploration = 2*bw^2/pi;
fprintf('There are done %i explorations convering %f sr in each exploration\n', ...
        explorations_needed, sr_covered_exploration);
    
% Exploration time
t_exploration = T/explorations_needed;
speed_signal = c;     % (m/us)
t_minimum = 2*R(end)/speed_signal;
fprintf('There are needed %f s per exploration. Each exploration takes %f s\n', ...
        t_minimum, t_exploration);
    
N_bits = getPhaseShifterBits(Nel_size, Nel_size, ...
    explorations_elevation, explorations_azimuth, elevation_angle, ...
    azimuth_angle);
fprintf('They are needed phase shifters of %i bits\n', N_bits);
    
% % Phase shifter bits
% N_bits = ceil(log2(explorations_needed));
% fprintf('For %i explorations, they are needed phase shifters of %i bits\n', ...
%         explorations_needed, N_bits);
    
    %% Question 2.1 Phase shifter bits
% In order to calculate the bits of the phase shifter, it is necessary to
% simulate the array with phase shifters of different bits. It is started
% with a fixed amount of bits (2), it is checked whether the beam formed is
% not less than 3 dB in the pointing direction than the ideal beam. If it
% failes, it is increased the number of bits until it is reached a number
% in which the theoretical and practical beams are similar.

% Step 1

%% Question 3
% calculate the bandwidth of the signal
tau = dc*t_exploration;
bw = 1/tau;
fprintf('BW = %f Hz\n', bw);

% Minimum noise figure of the lna
F_target = 3.8;             f_target = 10^(F_target/10);

F_circ = 0.2;               f_circ = 10^(F_circ/10);

F_lim = 0.8;                f_lim = 10^(F_lim/10);

f_max_lna = (f_target - f_circ - (f_lim - 1)*f_circ)/(f_circ*f_lim) + 1;
F_max_lna = 10*log10(f_max_lna);
fprintf('F max LNA = %f dB\n', F_max_lna);

% Complete noise figure
F_lna = 0.9;
G_lna = 14.8;

F_bpf = 2.2;

F_gr2106 = 1.1;
G_gr2106 = 20;

F_gr2374 = 1.8;
G_gr2374 = 14;

F_swi = 0.35;

F_ps = 13;

% They are put all the F and G together: the gain of the attenuators is its
% F with - sign. 
F_all = [F_circ F_lim F_lna F_bpf F_gr2106 F_gr2374 F_bpf F_swi F_ps];
G_all = [-F_circ -F_lim G_lna -F_bpf G_gr2106 G_gr2374 -F_bpf F_swi -F_ps];

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










