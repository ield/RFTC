function [J] = updateJ(J, Zin)
% Function created to update the J values when using the ADS
J([1 end]) = J([1 end])/sqrt(Zin);
end

