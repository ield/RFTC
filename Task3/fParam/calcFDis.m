function [F] = calcFDis(N, w, params, Z0, l, beta, alpha)
% Calculates the FParameters of each of the stages of the filter for each
% frequency when there are inverters and transmission lines
if(nargin == 6) 
    alpha = 0;
end

F = zeros(N+N+1, length(w), 2, 2);
for ii = 1:length(w)
    for jj = 1:2:N+N+1  % Inverters
        F(jj, ii, :, :) = [0 1i/params(jj); 1i*params(jj) 0];
    end
    for jj = 2:2:N+N+1  % Distributed element. See Pozar pag 190
        gamma = alpha + 1i*beta(ii);
        F(jj, ii, :, :) = [cosh(gamma*l) Z0*sinh(gamma*l); 1/Z0*sinh(gamma*l) cosh(gamma*l)];
    end
end
end

