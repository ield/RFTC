function [] = printCoupled(section, Z0e, Z0o, l, Z0, file)
%% Print title of the section
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, section);
fprintf(file, '\n');
fprintf(file, '\n');

%% Print l
fprintf(file, 'Z0e\n');
for ii = 1:length(Z0e)
    fprintf(file, '%i = %f \omega     ', ii, Z0e(ii));
end
fprintf(file, '\n');

fprintf(file, 'Z0o\n');
for ii = 1:length(Z0o)
    fprintf(file, '%i = %f \omega     ', ii, Z0o(ii));
end
fprintf(file, '\n');

fprintf(file, 'l = %d mm\n', l*1e3);

end

