function [] = printgk(g, file, text)
% Prints the gk elements of the low pass filter in file f
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, text);
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, 'Lumped elements\n');
C = g(1:2:end-1);   % Capacitors of the lumped elements
L = g(2:2:end-1);   % Inductances of the lumped elements.
% The last element is the matching resistance.

for ii = 1:length(C)
    fprintf(file, 'C%i = %d F\n', ii, C(ii));   
    fprintf(file, 'L%i = %d H\n', ii, L(ii));
end
fprintf(file, 'Rend = %d Ohm\n', 1/g(end));
end

