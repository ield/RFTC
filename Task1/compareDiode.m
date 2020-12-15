% Engineer: ield
% Comparisson of the beam lead pin diodes of ASI semiconductors to find the
% better match for a shunt configuration to implement in a phase shifter.
% Source: http://adsemi.com/pdf/diodes/BeamLeadPinDiodes.pdf
clear; close all;

Rs = 6.5; % Ohm
Cj = 0.02*1e-12;   % F

f0 = 10.25e9;               % Hz
f = [9:0.01:11] * 1e9;      % all frequency range (Hz)
Z0 = 50;                    % Ohm

% Calculate Zf and Zr. Approximations in Pozar pag 531 with L = 0
Zf = Rs*ones(size(f));
Zr = 1./(11*2*pi*f*Cj);

IL_on = -20*log10(abs(2*Zr./(2*Zr + Z0)));
IL_off = -20*log10(abs(2*Zf./(2*Zf + Z0)));


figure('Color',[1 1 1]);
set(gca, 'fontname','times new roman')

plot(f/1e9, IL_on); hold on;
plot(f/1e9, IL_off); hold on;

xlabel('Frequency (GHz)');
ylabel('IL (dB)');
legend('On', 'Off');

path = '../Task1_PhaseShifter/Report/Images/';
saveas(gca, [path, 'IL_diode'],'epsc');



