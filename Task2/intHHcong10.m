function [result] = intHHcong10(a, b, Emax)
% H*(H*) in z = 0 for m = 1, n = 0 = A^2bintegral(a)(cos^2(mpix/a))
eta = 120*pi;
% A. Knowing H = E/eta
A = Emax / eta;

% Solving the integral
intSize = 10;                         % Number of points in the integration
x = (-intSize/2:intSize/2)/intSize*a;   % Integration points
arg = (cos(pi*x/a)).^2;                  % Argument inside the integral
integ = sum(arg)/intSize;                     % Result of the integral

result = A^2*b*integ;

end

