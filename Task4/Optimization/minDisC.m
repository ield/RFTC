function dis = minDisC(C, f, Y_im)
% Calculates the distance between the curve of the admittance C and Y in
% the frequencies f
C = C*1e-12;
Yc_im = 2*pi*f*C;
dis = Yc_im-Y_im;
end

