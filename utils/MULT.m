function [SUM, DFA, DFE] = MULT(SNALP, SNBET, SALP, SBET,          ...
                                WAVLN, HA, XT, YT, ZT, ZN1, SIGMH, ...
                                XSPEC, XDIFF, ANTMOD)
% MULT  Compute monopulse sum and difference channel responses
%       including multipath effects over a spherical Earth.
% Subroutine MULT described in reference [1, p.40, 41]
%   Inputs:
%       SNALP, SNBET : Signal (measurement) azimuth and elevation angles [rad]
%       SALP,  SBET  : Steering (beam pointing) azimuth and elevation angles [rad]
%       WAVLN        : Wavelength [m]
%       HA           : Antenna height above the Earth surface [m]
%       XT, YT, ZT   : Cartesian coordinates of the target [m]
%       ZN1          : Relative dielectric constant of the reflecting surface
%       SIGMH        : RMS surface height (surface roughness) [m]
%       XSPEC        : Specular reflection weighting coefficient
%       XDIFF        : Diffuse reflection weighting coefficient
%       ANTMOD       : Elevation difference channel beam Configuration (antenna model)
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

    if nargin < 14, ANTMOD = 1; end
    AE = 8.5e6; % Effective Earth radius [m]
    R = sqrt(XT * XT + YT * YT + ZT * ZT);
    XL = sqrt(R * R - ZT * ZT);
%     XLT = sqrt(XT * XT + YT * YT);
    
    %% Determine location of multipath point in target co-ordinated
    Z = ZT + AE + HA;
    HT = sqrt(XL * XL + Z * Z) - AE;
    S = AE * acos(Z / sqrt(XL * XL + Z * Z));   % Arc length between antenna and target
    P = 2.e0 * sqrt(AE / 3.e0 * (HT + HA) + ((S / 2.e0)^2) / 3.e0);
    PHI = acos(2.e0 * AE * S * (HA - HT) / P^3);
    
    R1 = S / 2.e0 + P * cos((PHI + pi) / 3.e0); % Arc length from antenna to reflection point
    R2 = S - R1;                                % Arc length from reflection point to target
    ZL = AE * cos(R1 / AE); 
    XLM = sqrt(AE * AE - ZL * ZL);
    
    % Reflection point coordinates (XM, YM, ZM)
    ZM = ZL - AE - HA; 
    XM = XLM * cos(atan2(YT, XT));
    YM = XLM * sin(atan2(YT, XT));
    
    if AE * AE / ZL >= (AE + HA)
       fprintf('%s\n', 'Warning: Check the reflection point.');
    end
    
    %% Convert to Radar co-ordinates
    [RID1, SMULP, SMUET] = XYZ2RAzEl(XM, YM, ZM);

    %% Determine path length difference and critical angle
    RID = RID1 + sqrt((XT - XM)^2 + (YT - YM)^2 + (ZT - ZM)^2);
    ALPM = 2.e0 * pi * (RID - R) / WAVLN; % alpha
    PSI2 = atan(XM / ZL) + acos((abs(XM - XT)) / sqrt((XM - XT)^2 + (ZM - ZT)^2)); % Grazing angle
%     PSI2 = atan(XLM / ZL) + acos((abs(XLM - XLT)) / sqrt((XLM - XLT)^2 + (ZM - ZT)^2)); % Grazing angle
    PSIC = asin(sqrt(1.e0 / (ZN1 + 1.e0))); % ocean and soil 10 deg. (desert 30 deg.)
    if PSI2 <= PSIC, ALPM = ALPM + pi; end
    
    %% Calculate reflection coeff at multipath point
    TERM1 = ZN1 * sin(PSI2);
    TERM2 = sqrt(ZN1 - (cos(PSI2)^2));
    REFC = abs((TERM1 - TERM2) / (TERM1 + TERM2));                        % Rho_v
    D = 1.e0 / sqrt(1.e0 + 4.e0 * R1 * R2 / (AE * S * sin(2.e0 * PSI2))); % Divergence factor
    DELPHI = 4.e0 * SIGMH * sin(PSI2) / WAVLN;                            % DeltaPhi
    REFCS = sqrt(exp(-(DELPHI * pi)^2));
    
    %% Pick Rayleigh distributed sample for diffuse
    if DELPHI > 0.25e0
        SIGMD = 0.35e0;
    else
        SIGMD = 0.35e0 / 0.25e0 * DELPHI;
    end
    SIGMD = SIGMD / sqrt(2.e0);
    
    V1V2 = SIGMD * randn(2, 1);          % Normal distribution
    REFCD = sqrt(V1V2(1)^2 + V1V2(2)^2); % Amplitude, Rayleigh distribution
    ALPMD = rand(1) * 2.e0 * pi;         % Phase, uniform distribution.
    
    XIDS = REFC * D * REFCS * (cos(ALPM) + 1j * sin(ALPM)); % Specular
    XIDD = REFCD * REFC * (cos(ALPMD) + 1j * sin(ALPMD));   % Diffuse
    
    %% Calculate direct and indirect path antenna gain
    [SUMD, DFAD, DFED] = ANTENA(SNALP, SNBET, SALP, SBET, WAVLN, ANTMOD);
    [SUMID, DFAID, DFEID] = ANTENA(SMULP, SMUET, SALP, SBET, WAVLN, ANTMOD);
    
    %% Combine direct and indirect signals
    SUM = SUMD + XSPEC * SUMID * XIDS + XDIFF * SUMID * XIDD;
    DFA = DFAD + XSPEC * DFAID * XIDS + XDIFF * DFAID * XIDD;
    DFE = DFED + XSPEC * DFEID * XIDS + XDIFF * DFEID * XIDD;
end

