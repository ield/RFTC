function [] = plotQ(w, s11, s21, Q, path, file)
% Plots a normalized pass function
figure('Color',[1 1 1]);
% Set position of the plot
x0=500;
y0=500;
width=1300;
height=300;
set(gcf,'position',[x0,y0,width,height])
%%
subplot(1, 3, 1)
Legend = cell(length(Q)+1, 1);

data = abs(s11(1, :));
plot(w, 20*log10(data));
Legend{1}=strcat('Q = \infty');
hold on;

for ii = 2:length(Q)+1
    data = abs(s11(ii, :));
    plot(w, 20*log10(data), '--');
    Legend{ii}=strcat('Q = ', num2str(Q(ii-1)));
    hold on;
end

legend(Legend, 'location', 'northeast');

xlabel('Frequency [GHz]');
ylabel('s11 (dB)');
%%
subplot(1, 3, 2)
Legend = cell(length(Q)+1, 1);

data = abs(s21(1, :));
plot(w, 20*log10(data));
Legend{1}=strcat('Q = \infty');
hold on;

for ii = 2:length(Q)+1
    data = abs(s21(ii, :));
    plot(w, 20*log10(data), '--');
    Legend{ii}=strcat('Q = ', num2str(Q(ii-1)));
    hold on;
end

legend(Legend, 'location', 'south');

xlabel('Frequency [GHz]');
ylabel('s21 (dB)');
%%
subplot(1, 3, 3);
Legend = cell(length(Q)+1, 1);

data = s21(1, :)
phase = unwrap(angle(data));
gDelay = -diff(phase)./(diff(w)*2*pi);
plot(w(2:end), gDelay);

Legend{1}=strcat('Q = \infty');
hold on;

for ii = 2:length(Q)+1
    phase = unwrap(angle(s21(ii, :)));
    gDelay = -diff(phase)./(diff(w)*2*pi);

    plot(w(2:end), gDelay, '--');
    Legend{ii}=strcat('Q = ', num2str(Q(ii-1)));
    hold on;
end

legend(Legend, 'location', 'northeast');

xlabel('Frequency [GHz]');
ylabel('Group Delay(ns)');
ylim([0 4]);

%%
saveas(gca, [path, file],'epsc');
hold off;

end
