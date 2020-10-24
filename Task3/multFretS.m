function [s11,s21] = multFretS(lenPar, lenW, F)

    % Step Multiply all the F param of a given w (because they are
    % connected in series
    Ftot = zeros(lenW, 2, 2);
    res = zeros(2, 2);  % Accumulated result
    resp = zeros(2, 2); % Renewed result
    for ii = 1:lenW
        res(1, 1) = F(1, ii, 1, 1);
        res(1, 2) = F(1, ii, 1, 2);
        res(2, 1) = F(1, ii, 2, 1);
        res(2, 2) = F(1, ii, 2, 2);

        for jj = 2:lenPar
            resp(1, 1) = F(jj, ii, 1, 1);
            resp(1, 2) = F(jj, ii, 1, 2);
            resp(2, 1) = F(jj, ii, 2, 1);
            resp(2, 2) = F(jj, ii, 2, 2);
            res = res * resp;
        end
        Ftot(ii, 1, 1) = res(1, 1);
        Ftot(ii, 1, 2) = res(1, 2);
        Ftot(ii, 2, 1) = res(2, 1);
        Ftot(ii, 2, 2) = res(2, 2);
    end
    % Now we have a matrix of 3D where
    %   D1 = w: as if there was a matrix for each w.
    %   D3, D4: F param ([ABCD]) of the combined elements

    % Step 3.3. Convert to S param
    s11 = zeros(1, lenW);
    s21 = zeros(1, lenW);
    for ii = 1:lenW
        A = Ftot(ii, 1, 1);
        B = Ftot(ii, 1, 2);
        C = Ftot(ii, 2, 1);
        D = Ftot(ii, 2, 2);

        s11(ii) = (A+B-C-D)/(A+B+C+D);
        s21(ii) = 2/(A+B+C+D);
    end



end

