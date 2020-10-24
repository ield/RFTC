function [] = plotQ(w, sPar, Q, path, file, xla, yla)
% Plots a normalized pass function
figure('Color',[1 1 1]);
Legend = cell(length(Q)+1, 1);

plot(w, 10*log10(sPar(1, :)));
Legend{1}=strcat('Q = \infty');
hold on;

for ii = 2:length(Q)+1
    plot(w, 10*log10(sPar(ii, :)), '--');
    Legend{ii}=strcat('Q = ', num2str(Q(ii-1)));
    hold on;
end

legend(Legend)

xlabel(xla);
ylabel(yla);

saveas(gca, [path, file],'epsc');
hold off;

end
