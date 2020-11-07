function [] = plotcvsl(x, y, ylab, path, fileName)
figure('Color',[1 1 1]);
plot(x, y, 'k', 'LineWidth', 1); hold on;
xlabel('Patch length (mm)');
ylabel(ylab);
saveas(gca, [path, fileName],'epsc');
hold off;
end

