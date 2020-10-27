function [] = plotNormalized(w, s11, s21, path, file, xla, yla)
% Plots a normalized pass function
figure('Color',[1 1 1]);
% Set position of the plot
x0=500;
y0=500;
width=900;
height=300;
set(gcf,'position',[x0,y0,width,height])

% Adjust w to Hz
w = w/(2*pi);
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

phase = unwrap(angle(s21));
gDelay = -diff(phase)./(diff(w)*2*pi);
plot(w(2:end), gDelay);
hold on;

xlabel(xla);
ylabel(yla);
ylim([0 3]);

%%
saveas(gca, [path, file],'epsc');
hold off;

end

