function [] = plotNormalized(w, s11, s21, path, file, xla)
% Plots a normalized pass function
figure('Color',[1 1 1]);
% Set position of the plot
x0=500;
y0=500;
width=500;
height=300;
set(gcf,'position',[x0,y0,width,height])

% Adjust w to Hz
w = w/(2*pi);
%%

plot(w, 20*log10(abs(s11)));
hold on;
plot(w, 20*log10(abs(s21)));
 
xlabel(xla);
ylabel('dB');
legend('|s_1_1|', '|s_2_1|', 'location', 'southeast');
%%
saveas(gca, [path, file],'epsc');
hold off;

end

