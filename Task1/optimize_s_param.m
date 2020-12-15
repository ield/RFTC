function [error] = optimize_s_param(x, write)
    
    s12 = x(1);
    s22 = x(2);
    s23 = x(3);
    
    S_matrix = [[0.5 s12 1/sqrt(2)];[s12 s22 s23];[1/sqrt(2) s23 0]];
    
    error_matrix = S_matrix * conj(S_matrix) - eye(3);
    
    error = reshape(error_matrix,[],1); % convert matrix to column vector

    if(write)
        fprintf('s12 = %f; s22 = %f; s23 = %f. Error = %f\n', s12, s22, s23, error);
        
        fprintf('\n S \n');
        disp(S_matrix);
        
        fprintf('\n Error Matrix \n');
        disp(error_matrix);
        
    end
end

