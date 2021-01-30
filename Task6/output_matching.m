%% Develop the matching network
% The matching network is developed minimizing the rho of the matching, as
% done in the subject RFOT of the MSTC. The network can be formed of L or
% C, and depending on the side it is chosen either one or the other.
C_dcblock = 15e-6;                   % Capacitance of the dc block in F
% Zs = 50 + 1/(1j*win*C_dcblock);      % Input impedance of the source
Zs = 269.37*exp(1j*(4.054)*pi/180); % Input impedance in the output

% Zl = 71.74*exp(1j*(-90.233)*pi/180); % Input impedance of the transistor
% Zl = 37.739*exp(1j*(-43.569)*pi/180); % Input impedance of the transistor
Zl = 50;                            % Load impedance

match_elements = 2;     % Number of matching elements
xini = rand(1,match_elements)-0.5;
        

% The capacitances and indictances are expected to be in the range
%   1pf<c<1nf
%   1ph<l<1uh
% Therefore, it is seen that there is a gap between the minimum absolute
% value of l and c when computed the 'parameter p'. Therefore, that values
% are corrected as explained below.
% The limits are initially obtained for the capacitance and inductances values
% available in RS (https://es.rs-online.com/web/) and then corrected seeing
% which values lead to errors: l>2nH and c>10pf, c<0.5 pf

c_max = 100e-3;       % Value in F
c_min = 0.5e-12;      % Value in F
% Minimum abs(P) that can be obtained when C = c_min. Approx = 1.38e-7
min_abs_p_c = 1/(win*c_max);
% Minimum P that can be obtained when C = c_max. Approx = -6496
min_par = -1/(win*c_min);

l_min = 100e-12;      % Value in H
l_max = 100e-3;         % Value in H
% Minimum abs(P) that can be obtained when L = l_max. Approx = 1.5394
min_abs_p_l = win*l_min;
% Maximum P that can be obtained when L = L_max. Approx = 6e6
max_par = win*l_max;

lb = min_par*ones(1, match_elements);
ub = max_par*ones(1, match_elements);

sa_t=1000;
sa_rt=0.85;
sa_nt=5;
sa_ns=20;
[xopt, fopt]=simann(@(x) Matching_network_objective(x, fin, Zs, Zl, fin, ...
                   min_abs_p_c, min_abs_p_l),xini,lb,ub,...
                   sa_t,sa_rt,sa_nt,sa_ns,true);
               
rho_opt = matchingnetwork2(xopt, fin, Zs, Zl, fin);

fprintf('Simulated annheling results\n');
for ii = 1:length(xopt)
   if(xopt(ii) >= 0)
       par = xopt(ii) + min_abs_p_l;
       L = par / win*1e9;
       fprintf('L%i = %f nH; ', ii, L);
   else
       par = xopt(ii) - min_abs_p_c;
       C = -1 / (win*par)*1e12;
       fprintf('C%i = %f pF; ', ii, C);
   end
end
fprintf('RL = %f dB\n', -20*log10(abs(rho_opt)));

 