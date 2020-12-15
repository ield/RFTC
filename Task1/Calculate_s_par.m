clear;
xini = 0.5*rand(1, 3);

options=optimoptions('lsqnonlin');
options.Display='iter';
% lb = zeros(size(xini));
[xopt]=lsqnonlin(@(x) optimize_s_param(x, false),xini,[],[],options);

optimize_s_param(xopt, true);

	


