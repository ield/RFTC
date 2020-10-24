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
plotNormalized(w, s11, s21, path, 'LowPass');
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
plotNormalized(w, s11, s21, path, 'LPInverters1');
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
plotNormalized(w, s11, s21, path, 'LPInverters2');

%% Lowpass to bandpass transformation
f0 = 5e9;       % Central Frequency (Hz)
bw = 250e6;     % Bandwidth (Hz)
Z0 = 1;
w = 2*pi*(4e9:10e6:6e9);
f1 = f0-bw/2;
f2 = f0+bw/2;
f0 = sqrt(f1*f2);

% Step 1: Calculate the band pass parameters
% See Pozar pag 414
delta = bw/f0;
w0 = 2*pi*f0;

L = Z0*delta/w0;
C = 1/(w0*delta*Z0);

Y = 1./(1i*w) + 1i*w;
SP = sqrt(C/L);

% Step 2: Update inverters and parameters
J3 = ones(1, N+1);
for ii = 2:length(J3)-1
    J3(ii) = J2(ii)/SP;
end
J3(1) = J2(1)/sqrt(SP*Z0);
J3(N+1) = J2(N+1)/sqrt(SP*Z0);

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
[s11, s21] = multFretS(length(params), length(w), F);

% Step 5: Plot
plotNormalized(w/(2*pi), s11, s21, path, 'BPF');


