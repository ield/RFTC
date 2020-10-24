clear;
%% Low pass prototype
N = 3;
path = '../../Task3/Images/';   % Path to save the files
% Step 1: Calculate A max

RL = 20;                        % Return loss = 10log(Pret/Pin).
% Since Pin = Pret + Ptrans, and Att = 10log(Pin/Ptrans)
Amax = 10*log10(1/(1-10^(-RL/10)));

% Step 2: Calculate gk (Diapo 22)
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

% Step 3: Calculate S param
w = 0:0.01:5;
% Step 3.1: Calculate F param
F = zeros(N, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N % Capacitors
        F(jj, ii, :, :) = [1 0; 1i*g(jj)*w(ii) 1];
    end
    
    for jj = 2:2:N % Inductances
        F(jj, ii, :, :) = [1 1i*g(jj)*w(ii); 0 1];
    end
end

% Step 3.2: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11, s21] = multFretS(N, length(w), F);

% Step 4: Plot
plotNormalized(w, s11, s21, path, 'LowPass', 'Normalized \omega [rad/s]');
%% Ladder to parallel transformation with unitary inverters
% All inverters are introduced unitary and all the elements are transformed
% into capacitors with g = 1
J1 = ones(1, N+1);

params = zeros(1, N+N+1);   %j1 g1 j2 g2 j3 g3 j4
params(1:2:end) = J1;
params(2:2:end) = g(1:N);

% Step 1: Calculate F param
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Capacitors
        F(jj, ii, :, :) = [1 0; 1i*params(jj)*w(ii) 1];
    end
end

% Step 2: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11, s21] = multFretS(length(params), length(w), F);

% Step 4: Plot
plotNormalized(w, s11, s21, path, 'LPInverters1', 'Normalized \omega [rad/s]');
%% Ladder to parallel transformation with inverters
% Change the value of the inverters to have unitary capacitors
J2 = ones(1, N+1);
for ii = 2:length(J1)-1
    J2(ii) = J1(ii)/sqrt(g(ii-1)*g(ii));
end
J2(1) = 1/sqrt(g(1));
J2(N+1) = 1/sqrt(g(N));

params = zeros(1, N+N+1);   %j1 1 j2 1 j3 1 j4
params(1:2:end) = J2;
params(2:2:end) = ones(1, N);

% Step 1: Calculate F param
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Capacitors
        F(jj, ii, :, :) = [1 0; 1i*params(jj)*w(ii) 1];
    end
end

% Step 2: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11, s21] = multFretS(length(params), length(w), F);

% Step 3: Plot
plotNormalized(w, s11, s21, path, 'LPInverters2', 'Normalized \omega [rad/s]');

%% Lowpass to bandpass transformation, unnormalizing in Z and f (lumped elements)
f0 = 5e9;       % Central Frequency (Hz)
bw = 250e6;     % Bandwidth (Hz)
w = 2*pi*(4e9:10e6:6e9);
f1 = f0-bw/2;
f2 = f0+bw/2;
f0 = sqrt(f1*f2);
Z0 = 50;

% Step 1: Calculate the band pass parameters
% See Pozar pag 414
delta = bw/f0;
w0 = 2*pi*f0;

L = Z0*delta/w0;
C = 1/(w0*delta*Z0);

Y = 1./(1i*w*L) + 1i*w*C;

% Step 2: Update inverters and parameters
J3 = ones(1, N+1);
J3(2:end-1) = J2(2:end-1)/Z0;
J3([1 end]) = J2([1 end])/sqrt(Z0);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J3;
% Y1-3 are not inserted because they already depend on frequency. They will
% be used directly in the F param

% Step 3: Calculate F param
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Resonators
        F(jj, ii, :, :) = [1 0; Y(ii) 1];
    end
end

% Step 4: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11Lump, s21Lump] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotNormalized(w/(2*pi*1e9), s11Lump, s21Lump, path, 'BPF_Lumped', 'Frequency [GHz]');
fprintf('L = %f pH and C = %f pF \n', L*1e12, C*1e12);

%% Bandpass transformation, normalizing C, L

% Step 1: Calculate the Slope Parameter: see diapo 29
SP = sqrt(C/L);

% Step 2: Update inverters and parameters
J4 = ones(1, N+1);
J4(2:end-1) = J3(2:end-1)/SP;
J4([1 end]) = J3([1 end])/sqrt(SP);

params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J4;
Ynor = Y/SP;
% Y1-3 are not inserted because they already depend on frequency. They will
% be used directly in the F param

% Step 3: Calculate F param
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Resonators
        F(jj, ii, :, :) = [1 0; Ynor(ii) 1];
    end
end

% Step 4: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11, s21] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotNormalized(w/(2*pi*1e9), s11, s21, path, 'BPF_Normalized', 'Frequency [GHz]');

%% Transformation to distributed elements.
% Step 1: Calculate the distributed Slope Parameter: see radio notes on
% surf apuntes/t2 pag 33
% The change is done changing the resonators by lambda/2 waveguides of
% impezance Z0
SP = sqrt(C/L);
SPdis = pi/2/Z0;

% Step 2: Update inverters and parameters
J5 = ones(1, N+1);
J5(2:end-1) = J4(2:end-1)*SPdis;
J5([1 end]) = J4([1 end])*sqrt(SPdis);

params = zeros(1, N+N+1);   %j1 L1 j2 L2 j3 L3 j4
params(1:2:end) = J5;

% Step 3: Calculate F param
% Definition of F param of a line see Pozar pag 190
c = 3e8;
lambda = c/f0;

e = 3.2;    % Dielectric constant
lambdag = lambda/sqrt(e);
beta = w/c*sqrt(e);     % Phase constant see Pozar pag 20
l = lambdag/2;    % Length of the distributed element: designed for the f0

F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Distributed element. See Pozar pag 190
        F(jj, ii, :, :) = [cos(beta(ii)*l) 1i*Z0*sin(beta(ii)*l); 1i/Z0*sin(beta(ii)*l) cos(beta(ii)*l)];
    end
end

% Step 4: Calculate S param
% Multiply the F parameters of the cascade elements and transform them to s
% parameters
[s11Dis, s21Dis] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotTwo(w/(2*pi*1e9), s11Lump, s21Lump, s11Dis, s21Dis, path, 'BPF_Distributed', 'Frequency [GHz]');
%% Quality Factor: Lumped elements

% Step 1: Define new parameters: G and Q. See diapo 29
Q = [50 100 200 500 1000 2000]; %Several Q are going to be tested and compared
G = w0*C./Q;
% All sparameters obtained in the comparison of Q are stored in s11 and s21
% to be compared with the original lumped elements (Q = infinity)
s11Par = [abs(s11Lump)];
s21Par = [abs(s21Lump)];

% Step 2: Update inverters and parameters
params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J3;       % J as in lumped elements (J3)

% Step 3: Calculate F param
for kk = 1:length(Q)    %Several Q are going to be tested and compared
    F = zeros(N+N+1, length(w), 2, 2);
    for ii = 1:length(w)
        for jj = 1:2:N+N+1  % Inverters
            F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
        end

        for jj = 2:2:N+N+1  % Resonators
            F(jj, ii, :, :) = [1 0; Y(ii)+G(kk) 1]; %Y from lumped elements and G from diapo 29.
        end
    end

    % Step 4: Calculate S param
    % Multiply the F parameters of the cascade elements and transform them to s
    % parameters
    [s11, s21] = multFretS(length(params), length(w), F);
    s11Par = [s11Par; abs(s11)];
    s21Par = [s21Par; abs(s21)];
end
% Step 5: Plot
plotQ(w/(2*pi*1e9), s11Par, Q, path, 'BPF_Lumped_Q_s11', 'Frequency [GHz]', 's_1_1 [dB]');
plotQ(w/(2*pi*1e9), s21Par, Q, path, 'BPF_Lumped_Q_s21', 'Frequency [GHz]', 's_2_1 [dB]');

%% Quality Factor: Distributed elements

% Step 1: Define new parameters: G and Q. See diapo 29
e = 3.2;    % Dielectric constant
lambdag = lambda/sqrt(e);
Q = [50 100 200 500 1000 2000]; %Several Q are going to be tested and compared
alpha = pi*lambdag./(Q*lambda^2);
% All sparameters obtained in the comparison of Q are stored in s11 and s21
% to be compared with the original lumped elements (Q = infinity)
s11Par = [abs(s11Dis)];
s21Par = [abs(s21Dis)];

% Step 2: Update inverters and parameters
params = zeros(1, N+N+1);   %j1 Y1 j2 Y2 j3 Y3 j4
params(1:2:end) = J5;       % J as in lumped elements (J3)

% Step 3: Calculate F param
% Definition of F param of a line see Pozar pag 190
c = 3e8;
lambda = c/f0;
beta = w/c*sqrt(e);     % Phase constant see Pozar pag 20
l = lambdag/2;    % Length of the distributed element: designed for the f0

for kk = 1:length(Q)    %Several Q are going to be tested and compared
    F = zeros(N+N+1, length(w), 2, 2);
    for ii = 1:length(w)
        for jj = 1:2:N+N+1  % Inverters
            F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
        end

        for jj = 2:2:N+N+1  % Distributed element. See Pozar pag 190
            gamma = alpha(kk) + 1i*beta(ii);
            F(jj, ii, :, :) = [cosh(gamma*l) Z0*sinh(gamma*l); 1/Z0*sinh(gamma*l) cosh(gamma*l)];
        end
    end

    % Step 4: Calculate S param
    % Multiply the F parameters of the cascade elements and transform them to s
    % parameters
    [s11, s21] = multFretS(length(params), length(w), F);
    s11Par = [s11Par; abs(s11)];
    s21Par = [s21Par; abs(s21)];
end
% Step 5: Plot
plotQ(w/(2*pi*1e9), s11Par, Q, path, 'BPF_Dist_Q_s11', 'Frequency [GHz]', 's_1_1 [dB]');
plotQ(w/(2*pi*1e9), s21Par, Q, path, 'BPF_Dist_Q_s21', 'Frequency [GHz]', 's_2_1 [dB]');
