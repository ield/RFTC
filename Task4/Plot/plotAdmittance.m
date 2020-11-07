function [] = plotAdmittance(f, y, ylab, path, fileName)
figure('Color',[1 1 1]);
plot(f, real(y), 'LineWidth', 1); hold on;
plot(f, imag(y), 'LineWidth', 1); hold on;
legend('Re(Y(f))', 'Im(Y(f))', 'Location', 'northwest');
xlabel('f (GHz)');
ylabel(ylab);
saveas(gca, [path, fileName],'epsc');
hold off;
end

