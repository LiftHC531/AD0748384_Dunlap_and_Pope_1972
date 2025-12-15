function dAngle = DfSumRatio2Angle(sum, dif, config)
% DfSumRatio2Angle Convert difference-to-sum ratio to angular error.
%   Inputs:
%       sum    : Complex sum channel response (normalized) [volt]
%       dif    : Complex azimuth difference channel response [volt]
%       config : Elevation difference channel beam Configuration (antenna model) 
%
%   Outputs:
%       dAngle : Angular error estimate [rad]

    if nargin < 3, config = 1; end
    
    switch config
        case 1
            k2 = 0.01208e0; % Monopulse slope (rad per unit ratio)
        otherwise
            k2 = 0.02360e0; % Monopulse slope (rad per unit ratio)
    %         k2 = 0.01630e0; % Monopulse slope (rad per unit ratio)
    end
    
    difSumRatio = (real(sum) * real(dif) + imag(sum) * imag(dif)) / ...
                  (real(sum) * real(sum) + imag(sum) * imag(sum));
    dAngle = k2 * difSumRatio; % Angular error estimate [rad]
end