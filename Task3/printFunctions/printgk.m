function [] = printgk(N, g, file)
% Prints the gk elements of the low pass filter in file f
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, 'Low pass prototype\n');
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, 'Lumped elements\n');
for ii = 1:N+1
    fprintf(file, 'g%i = %d\n', ii, g(ii));    
end
end

