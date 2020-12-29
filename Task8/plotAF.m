N = 10;
path = '../../Task8_ArrayFactor/Images/';

% Plot a general overview of the array factor function
set(0, 'DefaultAxesFontName', 'Times New Roman');
figure('Color',[1 1 1]);
fplot(@(x) abs(sin(N/2*x)/sin(x/2)), [-3*pi 3*pi], 'k');

ylabel('AF');
yticks([0 N])
yticklabels({'0','N'})

xlabel('\psi');
xticks([-2*pi 0 2*pi])
xticklabels({'-2\pi','0','2\pi'})

saveas(gca, [path, 'af'],'epsc');

%% Plot the array factor function with different limits
limits = [0.5 1 2];     % The \psi limits (maximum and minimum) for a given separation
separation = {'Np_{\lambda/4}'; 'Np_{\lambda/2}'; 'Np_{\lambda}'};
xmin =  {'-\pi/2'; '-\pi'; '-2\pi'};
xmax =  {'\pi/2'; '\pi'; '2\pi'};
name = {'lambda_4', 'lambda_2', 'lambda'};
nPlots = 3;             % Number of plots to compare

for ii = 1:nPlots
    
%     subplot(1, nPlots, ii)
    figure('Color',[1 1 1]);
%     set(gcf,'units','normalized','outerposition',[0 0 0.33 1])

    fplot(@(x) abs(sin(N/2*x)/sin(x/2)), limits(ii)*[-pi pi], 'k');

    ylabel('AF');
    yticks([0 N])
    yticklabels({'0',separation{ii}})

    xlabel('\psi');
    xticks(limits(ii)*[-pi 0 pi])
    xticklabels({xmin{ii},'0',xmax{ii}})

    % Insert the \theta axis
    ax1 = gca;
    ax1_pos = ax1.Position; 
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none');
    x2 = linspace(0, pi, 100);
    y2 = zeros(size(x2));
    line(x2,y2,'Parent',ax2,'Color','k')

    yticks([])
    ylim([2 2.5]);

    xlabel('\theta');
    xticks([0 pi/2 pi])
    xticklabels({'0','\pi/2', '\pi',})
    xlim([0 pi]);
    
    filename = strcat('af_d_', name{ii});
    saveas(gca, [path, filename],'epsc');

end

close all;
%% Attempt to solve the integral
N = 10;
m = linspace(0.01, 3, 500);

psi_min = -2*pi*m;
psi_max = 2*pi*m;

psi_res = 180;      % Number of samples in psi
AF = zeros(length(m), psi_res); % All array factor samples will be here
allIntegralResults = zeros(size(m));

for ii = 1:length(m)
    psi = linspace(psi_min(ii), psi_max(ii), psi_res);    %180 psi samples of af function
    af = abs(sin(N/2*psi)./sin(psi/2));
    AF(ii,:) = af;
    allIntegralResults(ii) = sum(af)/psi_res;
end

d_05_pos = find(m>=0.5, 1);
p = allIntegralResults(d_05_pos)./allIntegralResults;

figure('Color',[1 1 1]);
plot(m, p, 'k');
ylabel('p');
xlabel('m');


% saveas(gca, [path, 'p_m'],'epsc');

figure('Color',[1 1 1]);
plot(m, N*p, 'k');
ylabel('Maximum directivity');
yticks([0 N])
yticklabels({'0','N'})

xlabel('m');

saveas(gca, [path, 'max_dir'],'epsc');

% Polynomial approximation of the directivity
orderPol = 30;
pol = polyfit(m,p,orderPol)
polApproximation = zeros(size(p));

fprintf('\n');
for ii = 1:orderPol+1
    power = orderPol - ii+1;
    fprintf('%f m^%i + ', pol(ii), power);  % Prints the polynomial
    polApproximation = polApproximation + pol(ii)*m.^power;
end
figure('Color',[1 1 1]);
plot(m, N*p, 'b'); hold on;
plot(m, N*polApproximation, 'r--');

ylabel('p');
xlabel('m');
yticks([0 N])
yticklabels({'0','N'})

legend('Maximum Directivity', 'Polynomial approximation', 'Location', 'Best');

saveas(gca, [path, 'pol_approx'],'epsc');
