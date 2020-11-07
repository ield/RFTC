function dis = minDisC(C, f, Y_im)
% Calculates the distance between the curve of the admittance C and Y in
% the frequencies f
Yc_im = 2*pi*f*C;
dis = sum((Yc_im-Y_im).^2);
end

