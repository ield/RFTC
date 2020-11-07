function [] = plotAdmittanceAndC(f, y, yc, ylab, path, fileName)
figure('Color',[1 1 1]);
plot(f, real(y), 'LineWidth', 1); hold on;
plot(f, imag(y), 'LineWidth', 1); hold on;

plot(f, real(yc), 'LineWidth', 1); hold on;
plot(f, imag(yc), 'LineWidth', 1); hold on;

legend('Re(Y(f))', 'Im(Y(f))', 'Re(Y_c(f))', 'Im(Y_c(f))', 'Location', 'northwest');
xlabel('f (GHz)');
ylabel(ylab);
saveas(gca, [path, fileName],'epsc');
hold off;
end

