% Figure 11 as described in reference [1, p.15, 26]
%   Reference:
%       [1] S. O. Dunlap and B. E. Pope, "Digital Simulation of Monopulse 
%           Angle Tracking with Multipath Propagation," Tech. Rep. AD0748384, 
%           U.S. Army Missile Command, Advanced Sensors Directorate, 
%           Redstone Arsenal, AL, USA, May 31, 1972.
clc; clear; close all;
tic;

addpath("./utils");

ft = 'Times';     % Times New Roman
ft_size = 14;     % Font size
set(0, 'defaultAxesFontName', ft, 'defaultTextFontName', ft,      ...
    'defaultUipanelFontName', ft, 'defaultUicontrolFontName', ft, ...
    'defaultUitableFontName', ft);

imageName1 = './figures/DP72_Fig11.png';
%%
c = 3.e8;                                  % Speed of light [m s-1]
fc = 8.e9;                                 % X-band operational frequency [Hz]
lambda = c / fc;                           % Wavelength [m]
specCoeff = 0.95e0;                        % Specular reflection coefficient
diffCoeff = 0.e0;                          % Diffuse reflection coefficient
Hr = 3.e0;                                 % Antenna height [m]
antennaModel = 1;                          % Antenna model
relDielectricConst = [3.e0, 15.e0, 30.e0]; % Relative dielectric constant of the reflecting surface
sigmaH = 0.e0;                             % RMS surface height [m]

numSamples = 1000;
% Target state initialization
Tar.Vx = 141.e0 * ones(1, numSamples);
Tar.Vy = 141.e0 * ones(1, numSamples);
Tar.Vz = 0.e0 * ones(1, numSamples);
Tar.X = zeros(1, numSamples); Tar.Y = zeros(1, numSamples);
Tar.Z = zeros(1, numSamples); Tar.R = zeros(1, numSamples);
Tar.Az = zeros(1, numSamples); Tar.El = zeros(1, numSamples);
Tar.X(1) = 200.e0;
Tar.Y(1) = 200.e0;
Tar.Z(1) = 50.e0;
[Tar.R(1), Tar.Az(1), Tar.El(1)] = XYZ2RAzEl(Tar.X(1), Tar.Y(1), Tar.Z(1));

% Radar state variables
Rdr.R = zeros(1, numSamples); Rdr.Az = zeros(length(relDielectricConst), numSamples); 
Rdr.El = zeros(length(relDielectricConst), numSamples);
Rdr.measAz = zeros(1, numSamples); Rdr.measEl = zeros(1, numSamples);
Rdr.errEl = zeros(length(relDielectricConst), numSamples); 
Rdr.errEl(:, 1) = NaN;

dt = 0.1e0;   % Data rate 10 Hz
stdR = 0.e0;  % Range measurement noise standard deviation [m]
stdAz = 0.e0; % Azimuth measurement noise standard deviation [rad]
stdEl = 0.e0; % Elevation measurement noise standard deviation [rad]
for k = 1 : length(relDielectricConst)
    Rdr.R(1) = Tar.R(1);
    Rdr.Az(k, 1) = Tar.Az(1);
    Rdr.El(k, 1) = Tar.El(1);
    for i = 2 : numSamples
        if k == 1
            Tar.X(i) = Tar.X(i - 1) + Tar.Vx(i) * dt;
            Tar.Y(i) = Tar.Y(i - 1) + Tar.Vy(i) * dt;
            Tar.Z(i) = Tar.Z(i - 1) + Tar.Vz(i) * dt;
            [Tar.R(i), Tar.Az(i), Tar.El(i)] = XYZ2RAzEl(Tar.X(i), Tar.Y(i), Tar.Z(i));
            Rdr.R(i) = Tar.R(i) + stdR*randn(1);
            Rdr.measAz(i) = Tar.Az(i) + stdAz * randn(1);
            Rdr.measEl(i) = Tar.El(i) + stdEl * randn(1);
        end

        [Sum, DfAz, DfEl] = MULT(Rdr.measAz(i), Rdr.measEl(i), ...
                                 Rdr.Az(k, i - 1), Rdr.El(k, i - 1), lambda, Hr, ...
                                 Tar.X(i), Tar.Y(i), Tar.Z(i), ...
                                 relDielectricConst(k), sigmaH, ...
                                 specCoeff, diffCoeff, antennaModel);
        dAz = DfSumRatio2Angle(Sum, DfAz);
        dEl = DfSumRatio2Angle(Sum, DfEl);

        % Beam steering update
        Rdr.Az(k, i) = Rdr.Az(k, i - 1) + dAz;
        Rdr.El(k, i) = Rdr.El(k, i - 1) + dEl;
        
        % Compute tracking error
        Rdr.errEl(k, i) = 1.e3 * (Rdr.El(k, i) - Tar.El(i));
    end
end

Tar.ElDeg = 180.e0 / pi * Tar.El;

%% Plot
figure; set(gcf, 'Color', 'w');
plot(Tar.ElDeg, Rdr.errEl(1, :), 'k-', 'LineWidth', 1.4); hold on;
plot(Tar.ElDeg, Rdr.errEl(2, :), 'b--', 'LineWidth', 1.4);
plot(Tar.ElDeg, Rdr.errEl(3, :), 'r-.', 'LineWidth', 1.4); hold off;

set(gca, 'XDir', 'reverse', 'XScale', 'log'); grid on;
xlim([0.0, 10.0]);
ylim([-10.0, 10.0]);
xLabels = flip([10.0, 5.0, 4.0, 3.0, 2.0, 1.5, 1.2, 1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.25]);
xticks(xLabels);
legend(strcat("\epsilon_r = ", num2str(relDielectricConst(1))), ...
       strcat("\epsilon_r = ", num2str(relDielectricConst(2))), ...
       strcat("\epsilon_r = ", num2str(relDielectricConst(3))), 'Location', 'southeast');
xlabel('Elevation angle [degrees]', 'FontSize', ft_size);
ylabel('Elevation channel tracking error [mrad]', 'FontSize', ft_size); grid on;
title('Figure 11. Multipath error specular reflection from smooth Earth', 'FontSize', ft_size);

width = 700;
height = 500;
set(gcf, 'position', [300, 100, width, height]);

% output the image
% print(gcf, '-dpng', '-r600', imageName1);

toc;
return;