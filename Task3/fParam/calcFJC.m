function [F] = calcFJC(N, w, params)
% Calculates the FParameters of each of the stages of the filter for each
% frequency when there are inverters and paralel capacitors
% Step 1: Calculate F param
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Capacitors
        F(jj, ii, :, :) = [1 0; 1i*params(jj)*w(ii) 1];
    end
end
end

