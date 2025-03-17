% Junsu Jang (junsu.jang94@gmail.com)
% 2025/01/20
function x = convertADC2V(x)
    % Convert the ADC digitial values to voltage
    ADrange    = (2^16)-1;
    voltRange = 5;
    x = x * (voltRange) / ADrange;  % Data in V
end