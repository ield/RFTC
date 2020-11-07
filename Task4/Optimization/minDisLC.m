function dis = minDisLC(LC, f, Y_im)
% Calculates the distance between the curve of the admittance LC and Y in
% the frequencies f
L = LC(1);
C = LC(2);
Ylc_im = 2*pi*f*C./(1-(2*pi*f).^2*L*C);
dis = sum((Ylc_im-Y_im).^2);
end

