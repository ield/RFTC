function [F] = calcFBPF(N, w, params, Y, G)
% Calculates the FParameters of each of the stages of the filter for each
% frequency when there are inverters and paralel resonances L||C
if(nargin == 4)
    G = 0;
end
F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    
    for jj = 2:2:N+N+1  % Resonators
        F(jj, ii, :, :) = [1 0; Y(ii)+G 1]; % G from Quality factor
    end
end
end

