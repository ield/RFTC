clear;
% Script to simulate the step by step design of a bpf from a lpf. All the
% steps follow:
% Step 1: Calculate the filter parameters: these are the value of the
%   lumped elements, transmision line or inverters that form the filter.
% Step 2: Calculate the F parameters (ABCD matrix) of each one of the
%   elements.
% Step 3: Multiply the F matrix of the cascade elements and obtain the s11
%   and s21 from the result.
% Step 4: Plot the s11, s21 and group velocity (all the plots are saved in 
%   the latex file) and/or print the values of the elements in the 
%   results.txt document.
%% Initial variables
N = 3;
path = '../../Task3/Images/';   % Path to save the files
resFile = 'results.txt';
file = fopen(resFile, 'w');

% Step 1: Calculate A max
RL = 20;                        % Return loss = 10log(Pret/Pin).
% Since Pin = Pret + Ptrans, and Att = 10log(Pin/Ptrans)
Amax = 10*log10(1/(1-10^(-RL/10)));
%% Low pass prototype
% Step 1: Calculate gk (Diapo 22)
beta = log(coth(Amax/17.37));
gamma = sinh(beta/(2*N));

index = 1:N;                    % Array 1-N to create values ak and bk

ak = sin((2*index-1)*pi/(2*N));
bk = gamma^2 + (sin(index*pi/N)).^2;

g = zeros(1, N+1);
g(1) = 2*ak(1)/gamma;
for ii = 2:N
    g(ii) = 4*ak(ii-1)*ak(ii)/(bk(ii-1)*g(ii-1));
end
g(N+1) = 1;

% Step 2:
w = 0:0.01:(2*pi);  % Normalized frequency
F = calculateFLPFNor(N, w, g);

% Step 3
[s11, s21] = multFretS(N, length(w), F);

% Step 4
plotNormalized(w, s11, s21, path, 'LowPass', 'Normalized f [Hz]', 'Group delay (s)');
printgk(N, g, file);
%% Ladder to parallel transformation with unitary inverters
% All inverters are introduced unitary and all the elements are transformed
% into capacitors with g = 1
% Step 1
J1 = ones(1, N+1);

params = zeros(1, N+N+1);   %j1 g1 j2 g2 j3 g3 j4
params(1:2:end) = J1;
params(2:2:end) = g(1:N);

% Step 2
F = calcFJC(N, w, params);
% Step 3
[s11, s21] = multFretS(length(params), length(w), F);

% Step 4: Plot and print values
plotNormalized(w, s11, s21, path, 'LPInverters1', 'Normalized f [Hz]', 'Group delay (s)');
printInvertersG('Ladder to paralel. Unitary inverters', N, F, g, file);

%% Normalization of capacitors (g values)
% Change the value of the inverters to have unitary capacitors
% Step 1
J2 = ones(1, N+1);
for ii = 2:length(J1)-1
    J2(ii) = J1(ii)/sqrt(g(ii-1)*g(ii));
end
J2(1) = 1/sqrt(g(1));
J2(N+1) = 1/sqrt(g(N));

params = zeros(1, N+N+1);   %j1 1 j2 1 j3 1 j4
params(1:2:end) = J2;
params(2:2:end) = ones(1, N);

% Step 2
F = calcFJC(N, w, params);
% Step 3
[s11, s21] = multFretS(length(params), length(w), F);

% Step 4
plotNormalized(w, s11, s21, path, 'LPInverters2', 'Normalized f [Hz]', 'Group delay (s)');
printInvertersG('Ladder to paralel. Unitary capacitors', N, F, params(2:2:end), file);
%% Lumped elements MATLAB
% Lowpass to bandpass transformation, unnormalizing in Z and f
f0 = 5e9;       % Central Frequency (Hz)
bw = 250e6;     % Bandwidth (Hz)
w = 2*pi*(4e9:10e6:6e9);
f1 = f0-bw/2;
f2 = f0+bw/2;
f0 = sqrt(f1*f2);
Z0 = 50;
Zin = 1;        % Input and output impedance. Set to 50 when using the data for the ADS

% Step 0: Calculate the band pass parameters
% See Pozar pag 414
delta = bw/f0;
w0 = 2*pi*f0;

L = Z0*delta/w0;
C = 1/(w0*delta*Z0);

Y = 1./(1i*w*L) + 1i*w*C;

% Step 1
J3 = ones(1, N+1);
J3(2:end-1) = J2(2:end-1)/Z0;
J3([1 end]) = J2([1 end])/sqrt(Z0);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J3;
% Y1-3 are not inserted because they already depend on frequency. They will
% be used directly in the F param

% Step 2
F = calcFBPF(N, w, params, Y);

% Step 3
[s11Lump, s21Lump] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotNormalized(w/(1e9), s11Lump, s21Lump, path, 'BPF_Lumped', 'Frequency [GHz]', 'Group delay (ns)');
% printInvertersLC('BPF_LumpedElements', N, F, L, C, Zin, file);
%% Lumped elements ADS
% Lowpass to bandpass transformation, unnormalizing in Z and f
Zin = 50;
% Step 1: Normalize inverters
J3_nor = updateJ(J3, Zin);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J3_nor;
% Y1-3 are not inserted because they already depend on frequency. They will
% be used directly in the F param

% Step 2
F = calcFBPF(N, w, params, Y);

% Step 5: Plot
printInvertersLC('BPF_LumpedElements ADS', N, F, L, C, Zin, file);
%% Bandpass transformation, normalizing C, L MATLAB
% Step 1: Calculate the Slope Parameter: see diapo 29
SP = sqrt(C/L);

% Step 2: Update inverters and parameters
J4 = ones(1, N+1);
J4(2:end-1) = J3(2:end-1)/SP;
J4([1 end]) = J3([1 end])/sqrt(SP);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J4;
Ynor = Y/SP;
Lnor = L*SP;
Cnor = C/SP;   
% Y1-3 are not inserted because they already depend on frequency. They will
% be used directly in the F param

% Step 3: Calculate F param
F = calcFBPF(N, w, params, Ynor);
% Step 4: Calculate S param
[s11, s21] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotNormalized(w/1e9, s11, s21, path, 'BPF_Normalized', 'Frequency [GHz]', 'Group delay (ns)');
% printInvertersLC('BPF Normalized Resonators', N, F, Lnor, Cnor, 1, file);

%% Bandpass transformation, normalizing C, L ADS
% Step 1
Zin = 50;
J4_nor = updateJ(J4, Zin);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J4_nor;

% Step 2: Calculate F param
F = calcFBPF(N, w, params, Ynor);

% Step 5: Plot
printInvertersLC('BPF Normalized Resonators ADS', N, F, Lnor, Cnor, Zin, file);
%% Transformation to distributed elements.
% Step 0: Calculate the distributed Slope Parameter: see radio notes on
% surf apuntes/t2 pag 33
% The change is done changing the resonators by lambda/2 waveguides of
% impezance Z0
SP = sqrt(C/L);
SPdis = pi/2/Z0;

c = 3e8;
lambda = c/f0;

e = 3.2;    % Dielectric constant
lambdag = lambda/sqrt(e);
beta = w/c*sqrt(e);     % Phase constant see Pozar pag 20
l = lambdag/2;    % Length of the distributed element: designed for the f0

% Step 1
J5 = ones(1, N+1);
J5(2:end-1) = J4(2:end-1)*SPdis;
J5([1 end]) = J4([1 end])*sqrt(SPdis);

params = zeros(1, N+N+1);   %j1 L1 j2 L2 j3 L3 j4
params(1:2:end) = J5;

% Step 3
F = calcFDis(N, w, params, Z0, l, beta);
% Step 4
[s11Dis, s21Dis] = multFretS(length(params), length(w), F);

% Step 5
plotTwo(w/(2*pi*1e9), s11Lump, s21Lump, s11Dis, s21Dis, path, 'BPF_Distributed', 'Frequency [GHz]');
% printDistributed('BPF Distributed elements', N, F, l, 1, file);
%% Transformation to distributed elements ADS.
% Step 1
J5_nor = updateJ(J5, Zin);

params = zeros(1, N+N+1);   %j1 L1 j2 L2 j3 L3 j4
params(1:2:end) = J5_nor;

% Step 3
F = calcFDis(N, w, params, Z0, l, beta);

% Step 5
printDistributed('BPF Distributed elements ADS', N, F, l, 1, file);
%% Quality Factor: Lumped elements
% Step 0: Define new parameters: G and Q. See diapo 29
Q = [50 100 200 500 1000 2000]; %Several Q are going to be tested and compared
G = w0*C./Q;
% All sparameters obtained in the comparison of Q are stored in s11 and s21
% to be compared with the original lumped elements (Q = infinity)
s11Par = [s11Lump];
s21Par = [s21Lump];

% Step 1: Update inverters and parameters
params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J3;       % J as in lumped elements (J3)

for kk = 1:length(Q)    %Several Q are going to be tested and compared
    % Step 3
    F = calcFBPF(N, w, params, Y, G(kk));
    % Step 4
    [s11, s21] = multFretS(length(params), length(w), F);
    s11Par = [s11Par; s11];
    s21Par = [s21Par; s21];
end
% Step 5
plotQ(w/(2*pi*1e9), s11Par, s21Par, Q, path, 'BPF_Lumped_Q_s11');
%% Quality Factor: Distributed elements
% Step 0: Define new parameters: G and Q. See diapo 29
e = 3.2;    % Dielectric constant
lambdag = lambda/sqrt(e);
Q = [50 100 200 500 1000 2000]; %Several Q are going to be tested and compared
alpha = pi*lambdag./(Q*lambda^2);
% All sparameters obtained in the comparison of Q are stored in s11 and s21
% to be compared with the original lumped elements (Q = infinity)
s11Par = [s11Dis];
s21Par = [s21Dis];

% Step 2: Update inverters and parameters
params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J5;       % J as in lumped elements (J3)


for kk = 1:length(Q)    %Several Q are going to be tested and compared
    % Step 2
    F = calcFDis(N, w, params, Z0, l, beta, alpha(kk));
    % Step 3
    [s11, s21] = multFretS(length(params), length(w), F);
    s11Par = [s11Par; s11];
    s21Par = [s21Par; s21];
end
% Step 5
plotQ(w/(2*pi*1e9), s11Par, s21Par, Q, path, 'BPF_Dist_Q_s11');
fclose(file);