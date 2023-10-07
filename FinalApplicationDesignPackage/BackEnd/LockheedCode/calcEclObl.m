function eps = calcEclObl (julianDate)
%--------------------------------------------------------------------------
% function eps = calcEclObl (julianDate)
%
% Calculate the true obliquity of the ecliptic plane
%--------------------------------------------------------------------------
%
% Inputs:
%     julianDate - Julian date
%     
% Outputs:
%     eps - True obliquity of the ecliptic plane [deg]
%
% References:     
%     Astronomical Formulae for Calculators, Jean Meeus, 
%     Willmann-Bell, Inc., Richmond, Virginia, 1988. (chapter 21)
%
%
% Date         Author                   Notes
% ----         ------                   -----
% 10/12/2005   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

TORAD = pi / 180.0;

days_from_J2000 = julianDate - 2451545.0;
T = days_from_J2000 / 36525.0;
T2 = T * T;
T3 = T2 * T;

% Sun's mean longitude [rad]
L = 280.4665 + 36000.7698 * T;
L = mod(L,360) * TORAD;  

% Moon's mean longitude [rad]
LP = 218.3165 + 481267.88134236 * T;
LP = mod(LP,360) * TORAD;   
                      
% Moon's argument of latitude (mean distance of the Moon from its 
% ascending node) [rad]

OM = 125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3/450000;
OM = mod(OM,360) * TORAD; 

% Mean obliquity of the ecliptic [deg]
eclip = (23.439291111111 - 46.8150/3600 * T - 0.00059/3600 * T2 ...
         + 0.001813/3600 * T3);

% Nutation in obliquity [deg]
deps = (9.20 * cos(OM) + 0.57 * cos(2*L) + 0.1 * cos(2*LP) - ...
        0.09 * cos(2*OM)) / 3600.0;
    
% True obliquity of the ecliptic [deg]
eps = eclip + deps;


 


