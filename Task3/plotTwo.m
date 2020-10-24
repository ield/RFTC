function [] = plotTwo(w, s11Lum, s21Lum, s11Dis, s21Dis, path, file, xla)
% Plots a normalized pass function
figure('Color',[1 1 1]);
% Set position of the plot
x0=500;
y0=500;
width=1100;
height=300;
set(gcf,'position',[x0,y0,width,height])

%%
subplot(1, 2, 1);
plot(w, 10*log10(abs(s11Lum)), '--');
hold on;
plot(w, 10*log10(abs(s21Lum)), '--');
hold on;
plot(w, 10*log10(abs(s11Dis)));
hold on;
plot(w, 10*log10(abs(s21Dis)));
hold on;

 
xlabel(xla);
ylabel('dB');
legend('|s_1_1| Lumped', '|s_2_1| Lumped', '|s_1_1| Distributed', '|s_2_1| Distributed', 'location', 'southeast');

%%
subplot(1, 2, 2);

phase = angle(s21Lum);
gDelay = -diff(phase)./diff(w);
plot(w(2:end), gDelay, '--');
hold on;
phase = angle(s21Dis);
gDelay = -diff(phase)./diff(w);
plot(w(2:end), gDelay);
hold on;

 
xlabel(xla);
ylabel('dB');
ylim([0 20]);
legend('Lumped', 'Distributed', 'location', 'northeast');



%%
saveas(gca, [path, file],'epsc');
hold off;

end

