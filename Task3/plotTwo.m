function [] = plotTwo(w, s11Nor, s21Nor, s11Dis, s21Dis, path, file, xla)
% Plots a normalized pass function
figure('Color',[1 1 1]);

plot(w, 10*log10(abs(s11Nor)), '--');
hold on;
plot(w, 10*log10(abs(s21Nor)), '--');
hold on;
plot(w, 10*log10(abs(s11Dis)));
hold on;
plot(w, 10*log10(abs(s21Dis)));
hold on;

 
xlabel(xla);
ylabel('dB');
legend('|s_1_1| Lumped', '|s_2_1| Lumped', '|s_1_1| Distributed', '|s_2_1| Distributed', 'location', 'northeast');

saveas(gca, [path, file],'epsc');
hold off;

end

