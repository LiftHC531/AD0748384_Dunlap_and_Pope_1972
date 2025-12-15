function [SUM, DFA, DFE] = ANTENA(SNALP, SNBET, SALP, SBET, WAVLN, CONFIG)
% This subroutine determines the response in free space of sin(x)/x pattern.
% Subroutine ANTENA described in reference [1, p.41]
%   Inputs:
%       SNALP, SNBET : Signal (measurement) azimuth and elevation angles [rad]
%       SALP,  SBET  : Steering (beam pointing) azimuth and elevation angles [rad]
%       WAVLN        : Wavelength [m]
%       CONFIG       : Elevation difference channel beam Configuration (antenna model)
%
%   Outputs:
%       SUM : Complex sum channel response (normalized) [volt]
%       DFA : Complex azimuth difference channel response [volt]
%       DFE : Complex elevation difference channel response [volt]
%
%   Reference:
%       [1] S. O. Dunlap and B. E. Pope, "Digital Simulation of Monopulse 
%           Angle Tracking with Multipath Propagation," Tech. Rep. AD0748384, 
%           U.S. Army Missile Command, Advanced Sensors Directorate, 
%           Redstone Arsenal, AL, USA, May 31, 1972.

    if nargin < 6, CONFIG = 1; end

    UAHA = 1.e0; UBHA = 1.e0; % Azimuth weighting parameters
    UAHE = 1.e0; UBHE = 1.e0; % Elevation weighting parameters
    
    % Parameters of NN and S are based on [1, p.46]
    NN = round(pi / (2.e0 * acot(45.8e0))); % ~72, cot(pi/(2*NN)) = 45.8
    S = 1.8e0 * WAVLN / pi; % Element spacing for rectangular array, P2 = pi*s/lambda = 1.8
    % S = 0.5e0 * WAVLN;
    
    XN = NN; % Number of elements per row or column
    P1 = pi * S * XN * 1.e0 / WAVLN;
    P2 = P1 / XN;
    P3 = WAVLN / (2.e0 * XN * S);
    
    F0 = 1.e0 / XN; % Normalization factor
    F1 = sin(P1 * (SNALP - SALP + P3)) / sin(P2 * (SNALP - SALP + P3)) * F0;
    F2 = sin(P1 * (SNALP - SALP - P3)) / sin(P2 * (SNALP - SALP - P3)) * F0;
    F3 = sin(P1 * (SNBET - SBET + P3)) / sin(P2 * (SNBET - SBET + P3)) * F0;
    F4 = sin(P1 * (SNBET - SBET - P3)) / sin(P2 * (SNBET - SBET - P3)) * F0;

    SUM = complex(0.5e0 * (F1 + F2) * (F3 + F4) / 0.81069946e0, 0.e0);
    DFA = (F2 - UBHA * F1) * (F3 + F4) / (0.81069946e0) * 0.5e0;
    DFE = (F4 - UBHE * F3) * (F1 + F2) / (0.81069946e0) * 0.5e0;

    % Elevation difference channel beam configuration 2 in [1, p.39]
    if CONFIG > 1
        X1 = SNBET - SBET; % [rad]
        X2 = abs(X1);
        DX = X2 * 180.e0 / pi; % [rad] -> [deg]

        if DX >= 0.e0 && DX <= 3.e0
            DE = 0.707e0 * sin(180.e0 / 3.e0 * X1);
        elseif DX > 3.e0 && DX <= 4.e0
            DE = 0.063e0 * sin(-180.e0 * (X1 - 3.e0 * pi / 180.e0));
        elseif DX > 4.e0
            DE = 0.01e0 * sin(180.e0 * (X1 - 4.e0 * pi / 180.e0));
        end
        DFE = DE;
    end
    
    DFA = complex(UAHA * DFA, 0.e0);
    DFE = complex(UAHE * DFE, 0.e0);
end