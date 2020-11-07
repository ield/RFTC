%Engineer: ield

function [collectedData] = textToSignal(text_s11, text_s31)
%% General Information
% Reads the s11 and s31 files generated with cst and returns an array of
% structs that each element:
%   l: length of the patch
%   f: frequencies evaluated
%   s11: s11 measures
%   s31: s31 measures

% It reads the files
file_s11 = fopen(text_s11);
data_s11 = textscan(file_s11,'%s');
fclose(file_s11);
data_s11 = str2double(data_s11{1}(1:end));

file_s31 = fopen(text_s31);
data_s31 = textscan(file_s31,'%s');
fclose(file_s31);
data_s31 = str2double(data_s31{1}(1:end));

% All the elements of the file have been scanned. Now the arrays are
% separated.

f1 = 0.1;   % First frequency sample
fend = 15;  % Last frequency sample

initialLength = 8;  % Minimum length of the patch, and first in the data
finalLength = 15;   % Maximum length of the patch, and last in the data
totalMeasures = (finalLength - initialLength)*2 + 1;    %Measures with different lengths

% collectedData = zeros(1, totalMeasures);     % Array of structures with variables length, f, s.

for ii = 1:totalMeasures
    collectedData(ii).length = initialLength + (ii-1)/2;     % Length of the patch
    
    % First located the values which correspond to the real and imaginary
    % measures of this length, therefore, it looks for the first two 0.1
    % and 15 in the data
    firstElem = find(data_s11 == f1, 2);    
    lastElem = find(data_s11 == fend, 2)+1;

    % Then the real section is the one in between the first 0.1 and 15 and
    % the imaginary section between the second
    realSection_s11 = data_s11(firstElem(1):lastElem(1));
    imagSection_s11 = data_s11(firstElem(2):lastElem(2));
    
    realSection_s31 = data_s31(firstElem(1):lastElem(1));
    imagSection_s31 = data_s31(firstElem(2):lastElem(2));

    f = realSection_s11(1:2:end);   % f is common for real and imaginary
    real_s11 = realSection_s11(2:2:end);
    imag_s11 = imagSection_s11(2:2:end);
    
    real_s31 = realSection_s31(2:2:end);
    imag_s31 = imagSection_s31(2:2:end);
    
    s11 = real_s11 + 1i*imag_s11;
    s31 = real_s31 + 1i*imag_s31;
    
    collectedData(ii).f = f;
    collectedData(ii).s11 = s11;
    collectedData(ii).s31 = s31;

    data_s11 = data_s11(lastElem(2)+1:end);    % Delete the values used of data, so that next it is searched for the next distance
    data_s31 = data_s31(lastElem(2)+1:end);
end



end

