function [] = plotNormalized(w, s11, s21, path, file, xla)
% Plots a normalized pass function
figure('Color',[1 1 1]);
% Set position of the plot
x0=500;
y0=500;
width=900;
height=300;
set(gcf,'position',[x0,y0,width,height])

%%
subplot(1, 2, 1);

plot(w, 10*log10(abs(s11)));
hold on;
plot(w, 10*log10(abs(s21)));
 
xlabel(xla);
ylabel('dB');
legend('|s_1_1|', '|s_2_1|', 'location', 'southeast');

%%
subplot(1, 2, 2);
% Plots group delay

phase = angle(s21);
gDelay = -diff(phase)./diff(w);
plot(w(2:end), gDelay);
hold on;

xlabel(xla);
ylabel('Group delay (s)');
ylim([0 20]);

%%
saveas(gca, [path, file],'epsc');
hold off;

end

