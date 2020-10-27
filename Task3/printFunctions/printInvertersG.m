function [] = printInvertersG(section, N, F, g, file)
%% Print title of the section
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, section);
fprintf(file, '\n');
fprintf(file, '\n');

%% Print F parameters of each inverter
fprintf(file, 'Inverters\n');
fpar = ['A' 'B'; 'C' 'D'];
for ii = 1:N+1  % For each inverter
    fprintf(file, 'J%i\n', ii);
    
    for jj = 1:2
        for kk = 1:2
            fprintf(file, '%c = %f + %fj;     ', fpar(jj,kk), real(F(2*ii-1, 1, jj, kk)), imag(F(2*ii-1, 1, jj, kk)));
        end
        fprintf(file, '\n');
    end
    fprintf(file, '\n'); 
end

%% Print g
fprintf(file, 'Lumped elements\n');

for ii = 1:length(g)
    fprintf(file, 'g%i = %d\n', ii, g(ii));    
end
end

