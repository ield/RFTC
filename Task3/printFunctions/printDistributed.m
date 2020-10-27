function [] = printDistributed(section, N, F, l, Z0, file)
%% Print title of the section
fprintf(file, '\n');
fprintf(file, '\n');
fprintf(file, section);
fprintf(file, '\n');
fprintf(file, '\n');

%% Print F parameters of each inverter for Z0
fprintf(file, 'Inverters for Z0 = %i\n', Z0);
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

%% Print l
fprintf(file, 'Microstrip lines\n');

fprintf(file, 'l = %d mm\n', l*1e3);

end

