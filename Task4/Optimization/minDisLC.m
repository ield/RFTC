function dis = minDisLC(LC, f, Y_im)
% Calculates the distance between the curve of the admittance LC and Y in
% the frequencies f
L = LC(1)*1e-9; % It is important to normalize
C = LC(2)*1e-12;
Ylc_im = 2*pi*f*C./(1-(2*pi*f).^2*L*C);
dis = Ylc_im-Y_im;
end

