function [] = plotNormalized(w, s11, s21, path, file, xla)
% Plots a normalized pass function
figure('Color',[1 1 1]);

plot(w, 10*log10(abs(s11)));
hold on;
plot(w, 10*log10(abs(s21)));
 
xlabel(xla);
ylabel('dB');
legend('|s_1_1|', '|s_2_1|', 'location', 'southeast');

saveas(gca, [path, file],'epsc');
hold off;

end

