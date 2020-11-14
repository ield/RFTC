function [F] = calculateFLPF(N, w, g)
% Calculates the FParameters of each of the stages of the filter for each
% frequency
F = zeros(N, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N % Capacitors
        F(jj, ii, :, :) = [1 0; 1i*g(jj)*w(ii) 1];
    end
    
    for jj = 2:2:N % Inductances
        F(jj, ii, :, :) = [1 1i*g(jj)*w(ii); 0 1];
    end
end
end

