% The figure 15 and 16 described in reference [1, pp.15 and 30-31]
%   Reference:
%       [1] S. O. Dunlap and B. E. Pope, "Digital Simulation of Monopulse 
%           Angle Tracking with Multipath Propagation," Tech. Rep. AD0748384, 
%           U.S. Army Missile Command, Advanced Sensors Directorate, 
%           Redstone Arsenal, AL, USA, May 31, 1972.
clc; clear; close all;
tic;

addpath("./utils");

ft = 'Times';     % Time News Roman
ft_size = 12;     % Font size
set(0, 'defaultAxesFontName', ft, 'defaultTextFontName', ft,      ...
    'defaultUipanelFontName', ft, 'defaultUicontrolFontName', ft, ...
    'defaultUitableFontName', ft);

imageName1 = './figures/DP72_Fig15_16.png';
%% Compute antenna pattern
c = 3.e8;        % Speed of light [m s-1]
fc = 8.e9;      % X-band operationl frequency [Hz]
lambda = c / fc; % Wavelength [m]

numSamples = 3000;
% Measurement angles [rad]
Az = 0.e0;                                          % SNALP in [1]
El = pi / 180.e0 * linspace(-5.0, 5.0, numSamples); % SNBET in [1]
% Steering angles [rad]
AzS = 0.e0;                                         % SALP in [1]
ElS = 0.e0;                                         % SBET in [1]

Sum1 = complex(zeros(1, numSamples), 0.e0);
DfAz1 = complex(zeros(1, numSamples), 0.e0);
DfEl1 = complex(zeros(1, numSamples), 0.e0);
DfEl2 = complex(zeros(1, numSamples), 0.e0);
for i = 1 : numSamples
    [Sum1(i), DfAz1(i), DfEl1(i)] = ANTENA(Az, El(i), AzS, ElS, lambda);
    [~, ~, DfEl2(i)] = ANTENA(Az, El(i), AzS, ElS, lambda, 2); % Antenna 2
end

ElDeg = 180.e0 / pi * El;
% Convert the volt to normalization dB
maxVal = max(abs(Sum1));
Sum1_dB = 20.e0 * log10(abs(Sum1) ./ maxVal);
DfEl1_dB = 20.e0 * log10(abs(DfEl1) ./ maxVal);
DfEl2_dB = 20.e0 * log10(abs(DfEl2) ./ maxVal);
%% Plot
figure; set(gcf, 'Color', 'w');
subplot(1, 2, 1);
plot(ElDeg, Sum1_dB, 'k--', 'LineWidth', 1.4); hold on;
plot(ElDeg, DfEl1_dB, 'r-', 'LineWidth', 1.4);
plot([-5.0, 5.0], [-3.0, -3.0], 'b-.', 'LineWidth', 1.0); hold off;
ylim([-50.0, 0.0]);
xticks(-5 : 1 : 5);

legend('\Sigma', '\Delta_{EL}', 'Location', 'southeast');
text(-4.5, -1.5, 'BW_{3dB} = 1.6^o');
xlabel('Degrees', 'FontSize', ft_size);
ylabel('Normalized antenna response [dB]', 'FontSize', ft_size); grid on;
title('Figure 15. Antenna 1 patterns', 'FontSize', ft_size);

subplot(1, 2, 2);
plot(ElDeg, Sum1_dB, 'k--', 'LineWidth', 1.4); hold on;
plot(ElDeg, DfEl2_dB, 'r-', 'LineWidth', 1.4);
plot([-5.0, 5.0], [-3.0, -3.0], 'b-.', 'LineWidth', 1.0); hold off;
ylim([-50.0, 0.0]);
xticks(-5 : 1 : 5);

xlabel('Degrees', 'FontSize', ft_size);
ylabel('Normalized antenna response [dB]', 'FontSize', ft_size); grid on;
title('Figure 16. Antenna 2 patterns', 'FontSize', ft_size);

width = 700;
height = 500;
set(gcf, 'position', [300, 100, width, height]);

% Output the image
% print(gcf, '-dpng', '-r600', imageName1);

toc;
return;