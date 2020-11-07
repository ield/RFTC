function [] = plotAdmittance(x, y, ylab, path, fileName)
figure('Color',[1 1 1]);
plot(x, real(y), 'LineWidth', 1); hold on;
plot(x, imag(y), 'LineWidth', 1); hold on;
legend('Re(Y(f))', 'Im(Y(f))', 'Location', 'northwest');
xlabel('f (GHz)');
ylabel(ylab);
saveas(gca, [path, fileName],'epsc');
hold off;
end

