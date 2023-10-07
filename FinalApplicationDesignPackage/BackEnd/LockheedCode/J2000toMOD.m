function tmatrix = J2000toMOD(julianDate)
%--------------------------------------------------------------------------
% function tmatrix = J2000toMOD(julianDate)
%
% Calculate the transformation matrix from J2000 to MOD frame
%--------------------------------------------------------------------------
%
% Inputs:
%     julianDate - Julian date
%     
% Outputs:
%     tmatrix - transformation matrix from J2000 to MOD frame
%
% References:     
%     Satellite Orbits Models Methods Applications,
%     Oliver Montenbruck and Eberhard Hill (chapter 5)
%     Astronomical almanac 1994 (page B18)
%
% Date         Author                   Notes
% ----         ------                   -----
% 06/12/2010   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

TORAD = pi / 180.0;

days_from_J2000 = julianDate - 2451545.0;
T = days_from_J2000 / 36525.0;
T2 = T * T;
T3 = T2 * T;

eta = (2306.2181/3600 * T + 0.30188/3600 * T2 + 0.017998/3600 * T3) * TORAD;
nu =  (2004.3109/3600 * T + 0.42665/3600 * T2 - 0.041833/3600 * T3) * TORAD; 
zeta = eta + (0.79280/3600 * T2 + 0.000205/3600 * T3) * TORAD;

tmatrix(1,1) = -sin(zeta) * sin(eta) + cos(zeta) * cos(nu) * cos(eta);
tmatrix(2,1) =  cos(zeta) * sin(eta) + sin(zeta) * cos(nu) * cos(eta);
tmatrix(3,1) =  sin(nu) * cos(eta);

tmatrix(1,2) = -sin(zeta) * cos(eta) - cos(zeta) * cos(nu) * sin(eta);
tmatrix(2,2) =  cos(zeta) * cos(eta) - sin(zeta) * cos(nu) * sin(eta);
tmatrix(3,2) = -sin(nu) * sin(eta);

tmatrix(1,3) = -cos(zeta) * sin(nu);
tmatrix(2,3) = -sin(zeta) * sin(nu);
tmatrix(3,3) =  cos(nu);



 

