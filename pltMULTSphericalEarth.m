% Schematic of multipath propagation over a spherical Earth
%   Reference:
%       [1] S. O. Dunlap and B. E. Pope, "Digital Simulation of Monopulse 
%           Angle Tracking with Multipath Propagation," Tech. Rep. AD0748384, 
%           U.S. Army Missile Command, Advanced Sensors Directorate, 
%           Redstone Arsenal, AL, USA, May 31, 1972.
%
%   Description:
%       This script reproduces the schematic of the spherical-Earth 
%       multipath geometry described in [1, p.5, 40]. It computes the 
%       antenna, target, and reflection-point locations, and plots the 
%       propagation paths including the direct and ground-reflected rays.
clc; clear; close all;
tic;

ft = 'Times';     % Time News Roman
ft_size = 14;     % Font size
set(0, 'defaultAxesFontName', ft, 'defaultTextFontName', ft,      ...
    'defaultUipanelFontName', ft, 'defaultUicontrolFontName', ft, ...
    'defaultUitableFontName', ft);

imageName1 = './figures/multipathModel_sphericalEarthSchematic.png';
%% Configuration
AE = 8.5e6; % Effective Earth radius [m]
numSamples = 1000;
theta = pi / 180.e0 * linspace(60.e0, 92.e0, numSamples);

HA = 10.e4; % Antenna height above Earth surface [m]
offset = AE + HA;
% Define the Earth surface profile
XE = AE * cos(theta);
ZE = AE * sin(theta) - offset;

% Antenna coordinates (XA, YA, ZA)
XA = 0.e0; YA = 0.e0;
ZA = 0.e0;

% Target coordinates (XT, YT, ZT)
XT = 30.e5; YT = 0.e0;
ZT = 10.e4;
R = sqrt(XT * XT + YT * YT + ZT * ZT);
XL = sqrt(R * R - ZT * ZT);

%% Compute spherical-Earth multipath model from [1, p.40]
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

%% Plot geometry
figure; set(gcf, 'Color', 'w');
plot(XE, ZE, 'k-', 'LineWidth', 1.4); hold on;
plot([XT, XA], [ZT, ZA], 'b-', 'LineWidth', 1.4); % Direct path
plot([XA, XM], [ZA, ZM], 'r-', 'LineWidth', 1.4); % Incident path
plot([XM, XT], [ZM, ZT], 'r-', 'LineWidth', 1.4); % Reflected path

% Auxiliary lines to the Earth's center
plot([XA, 0.0], [ZA, -offset], 'k--', 'LineWidth', 1.4);
plot([XT, 0.0], [ZT, -offset], 'k--', 'LineWidth', 1.4);
plot([XM, 0.0], [ZM, -offset], 'k--', 'LineWidth', 1.4);

plot([XA, XA], [ZA, Z - offset], '-.', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.4);
plot([XA, XL], [Z - offset, ZT], '-.', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.4);
plot([XA, XLM], [ZM, ZM], '-.', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1.4); hold off;

text(1.e6, 14.e4, 'X_L');
text(0.2e6, -1.6e5, 'X_{LM}');
text(-12.e4, 7.e4, 'Z_T');
text(-12.e4, -4.e4, 'H_A');
text(35.e4, -5.e5, 'A_E');
text(2.9e6, -3.e5, 'H_T');

xlim([-2.e5, 3.2e6]);
ylim([-9.5e5, 3.e5]);
xlabel('X [m]', 'FontSize', ft_size);
ylabel('Z [m]', 'FontSize', ft_size); grid on;
title('Schematic of multipath propagation over a spherical Earth', 'FontSize', ft_size);

width = 700;
height = 500;
set(gcf, 'position', [300, 100, width, height]);

% output the image
% print(gcf, '-dpng', '-r600', imageName1);

toc;
return;