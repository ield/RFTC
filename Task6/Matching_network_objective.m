function [error] = Matching_network_objective(x, f, Z_S, Z_L, f0, ...
                    min_abs_p_c, min_abs_p_l)
% Author: ield
% Transforms the rho given by matchingnetwrok2 to the error of the problem.
% The error is calculated as in a lsq nonlin problem.
% The error considers reducing rho.
% Type: 1 if the optimization, 2 if the optimization is local (lsqnonlin)
%% Transforms the values of x:
% The capacitances and indictances are expected to be in the range
%   1pf<c<1nf
%   1ph<l<1uh
% Therefore, it is seen that there is a gap between the minimum absolute
% value of l and c when computed the 'parameter p'. Therefore, that values
% are corrected as explained below
% 
% min_abs_p_c = Minimum abs(P) that can be obtained when C = c_max. Approx = -0.0065
% min_abs_p_l = Minimum abs(P) that can be obtained when L = l_max. Approx = 0.01539
%
% The values that are positive are increased min_abs_p_l and the values
% that are negative are decreased min_abs_p_c

x(x>=0) = x(x>=0) + min_abs_p_l;
x(x<0) = x(x<0) - min_abs_p_c;
%% Calculating the error
% There errors is due to have pho > 0
rho = matchingnetwork2(x, f, Z_S, Z_L, f0);

error = abs(rho);
end

