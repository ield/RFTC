function [D] = max_Dir_separation(m, a, alpha)
 % Returns the maximum directivity (D) as a function of the separation
 % between elements (m), the distribution a, and the progressive phase
 % alpha. It is used eq 5.18 of Cardama's book, knowing that kd = 2pim

    N = length(a);
    D = zeros(size(m));

    for ii = 1:length(m)   
        numerator = sum(a)^2;

        den_term1 = sum(a.^2);

        den_term2 = 0;
        for jj = 1:N-2
            for kk = jj+1:N-1
                den_term2 = den_term2 + 2*a(jj)*a(kk)*...
                    sin(2*pi*m(ii)*(jj-kk))/(2*pi*m(ii)*(jj-kk))*...
                    cos(alpha*(jj-kk));
            end
        end

    denominator = den_term1 + den_term2;    
    D(ii) = numerator/denominator;   
    end

end

