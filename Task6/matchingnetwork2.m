function [rho] = matchingnetwork2(params,f,Z_S,Z_L,f0)
%MATCHINGNETWORK computes the reflection coefficient of a 
%                complex-impedance matching network with 5 elements.
%   [rho] = MATCHINGNETWORK5(params,f,Z_S,Z_L)
%
%   params(5) = vector of reactances at f0
%   f = frequency (GHz)     (vector, same size as Z_S and Z_L)
%   Z_S = impedance (ohms)
%   Z_L = impedance (ohms)
%
%   rho = reflection coeffcient at frequencies f

omega0=2*pi*f0;
omega=2*pi*f;

if params(1)>=0
    Zp1=1j*(omega/omega0)*params(1);
else
    Zp1=1j*(omega0./omega)*params(1);
end
Yp1=1./Zp1;
if params(2)>=0
    Zs1=1j*(omega/omega0)*params(2);
else
    Zs1=1j*(omega0./omega)*params(2);
end

Zv=1./(Yp1+1./(Zs1+Z_L));
rho=(Zv-conj(Z_S))./(Zv+conj(Z_S));

end
