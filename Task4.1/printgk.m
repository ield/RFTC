function [] = printgk(g, file, text, elements, units)
% Prints the gk elements of the low pass filter in file f
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, text);
fprintf(file, '\n');
fprintf(file, '\n');
% The last element is the matching resistance.

for ii = 1:length(g)
    fprintf(file, '%s_%i = %f %s\n', elements, ii, g(ii), units);
end

end

