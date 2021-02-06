function [N_bits] = getPhaseShifterBits(M, N,...
                    explorations_elevation, explorations_azimuth, ...
                    elevation_angle, azimuth_angle)
% Question 2.1 Phase shifter bits
% In order to calculate the bits of the phase shifter, it is necessary to
% simulate the array with phase shifters of different bits. It is started
% with a fixed amount of bits (2), it is checked whether the beam formed is
% not less than 3 dB in the pointing direction than the ideal beam. If it
% failes, it is increased the number of bits until it is reached a number
% in which the theoretical and practical beams are similar.

% Inputs
% M:                        number of elements in x dimension
% N:                        number of elements in y dimension
% explorations_needed:      number of explorations that must be achieved
% explorations_elevation:   number of explorations in elevation
% explorations_azimuth:     number of explorations in azimuth
% elevation_angle:          [rad]
% azimuth_angle:            [rad]

%% Step 1. It is determined the elevations and azimuth angles in which the
% antenna is going to point
phi = linspace(-elevation_angle/2, elevation_angle/2, ...
                        explorations_elevation)+pi/2;
theta = linspace(-azimuth_angle/2, azimuth_angle/2, ...
                        explorations_azimuth);

%% Step 2. For every theta_0, phi_0 calculate betax and betay
beta_x_theoretical = zeros(explorations_azimuth, explorations_elevation);
beta_y_theoretical = zeros(explorations_azimuth, explorations_elevation);

for ii = 1:explorations_azimuth
    for jj = 1:explorations_elevation
        theta_0 = theta(ii);
        phi_0 = phi(jj);
        beta_x_theoretical(ii, jj) = -pi*sin(theta_0)*cos(phi_0);
        beta_y_theoretical(ii, jj) = -pi*sin(theta_0)*sin(phi_0);
    end    
end

% Check this part
% [azimuth_positions, elevation_positions] = meshgrid(azimuth_positions, elevation_positions);
% 
% size(azimuth_positions)
% size(elevation_positions)
% size(beta_y_theoretical)
% 
% 
% surf(azimuth_positions, elevation_positions, beta_x_theoretical');


for N_bits = 1:100
% for N_bits = 2:2
    N_bits_possible = true; % This is true until one of the AF is < -3
    %% Step 3. For every number of bits calculate the possible angles
    N_angles = 2^N_bits;
    available_phase = linspace(0, 2*pi, N_angles);
    AF = zeros(explorations_azimuth, explorations_elevation);
    
    %% Step 4. Calculte the phase of each array
    % For each position in azimuth and elevation calculate the ideal phase
    % shift and the nearest available phase shift
    for ii = 1:explorations_azimuth
        for jj = 1:explorations_elevation
            beta_x = beta_x_theoretical(ii, jj);
            beta_y = beta_y_theoretical(ii, jj);
            theta_0 = theta(ii);
            phi_0 = phi(jj);
            
            element_shifts = zeros(M, N);
            
            for mm = 1:M
                for nn = 1:N
                    % Be careful!! The theoretical phaseshift must be
                    % modulus 2*pi
                    theoretical_phase_shift_mn = ...
                        mod(mm*beta_x + nn*beta_y, 2*pi);
                    [~, pos_phase_shift] = min(abs(available_phase - ...
                        theoretical_phase_shift_mn));
                    element_shifts(mm, nn) = available_phase(pos_phase_shift);

                end
            end
%             % Check if the phase shift is chosen properly 
%             [mm, nn] = meshgrid(1:M, 1:N);
%             surf(mm, nn, element_shifts);
%             return;

            %% Step 5: Calculate the array factor
            % Using the phase shifts recently calculated
            
            for mm = 1:M
                for nn = 1:N
                    AF (ii, jj) = AF (ii, jj) + ...
                        exp(1j*((mm-1)*pi*sin(theta_0)*cos(phi_0) + ...
                        (nn-1)*pi*sin(theta_0)*sin(phi_0) + ...
                        element_shifts(mm, nn)));
                end
            end
            
            %% Step 6: check if the array factor is < -3
            % It is also normalized the array factor
            AF(ii, jj) = 10*log10(abs(AF(ii, jj)/(M*N)));
            if(AF(ii, jj) < -3)
                N_bits_possible = false;
            end            
        end 
%         if(~N_bits_possible) 
%             break;
%         end
    end
    
    if(N_bits_possible)
        fprintf('The phased array is achievable with %i bits.\n', N_bits);
        % plot the af
%         figure('Color',[1 1 1]);
        subplot(1, 2, 2);
        [theta_axis, phi_axis] = meshgrid(theta, phi);
        pcolor(theta_axis*180/pi, phi_axis*180/pi, AF');
        xlabel('\theta (º)');
        ylabel('\phi (º)');
        title(strcat('AF reduction (dB) for N_{bits} = ', num2str(N_bits)));
        colormap winter
        colorbar
        
%         path = '../../Task9_PhasedArray/Images/';
%         saveas(gca, [path, 'AF_reduction'],'epsc');

        return;
    else
        fprintf('The phased array is not achievable with %i bits.\n', N_bits);
        
        % Plor the array factor
        figure('Color',[1 1 1]);
        set(gcf,'position',[100,100,1000,400]);
        subplot(1, 2, 1);
        [theta_axis, phi_axis] = meshgrid(theta, phi);
        pcolor(theta_axis*180/pi, phi_axis*180/pi, AF');
        xlabel('\theta (º)');
        ylabel('\phi (º)');
        title(strcat('AF reduction (dB) for N_{bits} = ', num2str(N_bits)));
        colormap winter
        colorbar

%         break;
    end
    
end

end

