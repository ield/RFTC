function [Zte] = obtainZte(kc,f)
    c = 3e8;
    eta = 120*pi;
    
    k = 2*pi*f/c;
    beta = sqrt(k.^2 - kc^2);
    
    Zte = k*eta./beta;
    
end

