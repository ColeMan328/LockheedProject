function eclip = calcEclMeanObl (julianDate)
%--------------------------------------------------------------------------
% function eclip = calcEclMeanObl (julianDate)
%
% Calculate the mean obliquity of the ecliptic plane
%--------------------------------------------------------------------------
%
% Inputs:
%     julianDate - Julian date
%     
% Outputs:
%     eclip - Mean obliquity of the ecliptic plane [deg]
%
% References:     
%     Astronomical Formulae for Calculators, Jean Meeus, 
%     Willmann-Bell, Inc., Richmond, Virginia, 1988. (chapter 21)
%
%
% Date         Author                   Notes
% ----         ------                   -----
% 06/11/2010   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

days_from_J2000 = julianDate - 2451545.0;
T = days_from_J2000 / 36525.0;
T2 = T * T;
T3 = T2 * T; 

% Mean obliquity of the ecliptic [deg]
eclip = (23.439291111111 - 46.8150/3600 * T - 0.00059/3600 * T2 ...
         + 0.001813/3600 * T3);


 


