clear;
close all;

N = 10;
a = ones(1, N)/N;     % Feeding of each element

path = '../../Task8_ArrayFactor/Images/';

% Plot a general overview of the array factor function
set(0, 'DefaultAxesFontName', 'Times New Roman');
figure('Color',[1 1 1]);

psi = linspace(-3*pi, 3*pi, 10000);

AF = zeros(size(psi));

for ii = 1:length(AF)
    for jj = 1:N
        AF(ii) = AF(ii) + a(jj)*exp(1j*jj*psi(ii));
    end
end

plot(psi, abs(AF).^2, 'k');

ylabel('|AF|^2');

xlabel('\psi');
xticks([-2*pi 0 2*pi])
xticklabels({'-2\pi','0','2\pi'})

% saveas(gca, [path, 'af'],'epsc');

%% Different distributions with alpha = 1
% It is used eq 5.18 of Cardama's book, knowing that kd = 2pim
clear;
close all;

path = '../../Task8_ArrayFactor/Images/';

% Plot a general overview of the array factor function
set(0, 'DefaultAxesFontName', 'Times New Roman');
figure('Color',[1 1 1]);

N = 10;
m = linspace(0, 3, 1000);
alpha = 0;              % Phase progressive

% Uniform distribution
a = ones(1, N)/N;       % Feeding of each element
D_uniform = max_Dir_separation(m, a, alpha);

plot(m, D_uniform); hold on;

% Triangular distribution
a = zeros(1, N);
for ii = 0:N-1
    a(ii+1) = 1-abs((-(N-1)/2+ii)./(N/2));
end
a = a/sum(a);
D_triangular = max_Dir_separation(m, a, alpha);

plot(m, D_triangular); hold on;

% Cosine distribution
a = zeros(1, N);
H = 0.5;
for ii = 0:N-1
    elem = ii-(N-1)/2;
    a(ii+1) = 1+H*(cos(pi*elem/(N-1)))^2;
end
a = a/sum(a);
D_cosine = max_Dir_separation(m, a, alpha);

plot(m, D_cosine); hold on;

% Binomial distribution
a = zeros(1, N);
for ii = 0:N-1
    a(ii+1) = nchoosek(N-1,ii);
end
a = a/sum(a);
D_binomial = max_Dir_separation(m, a, alpha);

plot(m, D_binomial); hold on;

xlabel('m = d/\lambda');
ylabel('D');
legend('Uniform', 'Triangular', 'Cosine H = 0.5', 'Binomial');

% saveas(gca, [path, 'distributions_alpha0'],'epsc');

%% Uniform distribution with different alpha
% Progressive phase uniform distribution
clear;
close all;

path = '../../Task8_ArrayFactor/Images/';

% Plot a general overview of the array factor function
set(0, 'DefaultAxesFontName', 'Times New Roman');
figure('Color',[1 1 1]);

alpha = 0:15:90;
txt = cell(length(alpha),1);        % Legend of the graph

N = 10;
a = ones(1, N)/N;       % Feeding of each element
m = linspace(0, 1.5, 1000);

for ii = 1:length(alpha)
    D_phase = max_Dir_separation(m, a, alpha(ii)*pi/180);
    plot(m, D_phase); hold on;
    txt{ii}= sprintf('%i º',alpha(ii));
end
legend(txt);
xlabel('m = d/\lambda');
ylabel('D');

saveas(gca, [path, 'distributions_progrssive_phase'],'epsc');

