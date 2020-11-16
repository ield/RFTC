function [] = plotcvsl(N, Copts, lengths, Cres, Ldes, path, fileName)
% Copts and Cres must come in pF
figure('Color',[1 1 1]);

% Plots the lvsc, the inverse of the graph obtained in the previous lab
plot(Copts, lengths, 'k', 'LineWidth', 1); hold on;

% Plots the L obtained using interpolation
plot(Cres, Ldes, 'ro'); hold on;

xlabel('C (pF)');
ylabel('Patch Length (mm)');

fileName = [fileName, '_n', num2str(N)];
saveas(gca, [path, fileName],'epsc');
hold off;
end

