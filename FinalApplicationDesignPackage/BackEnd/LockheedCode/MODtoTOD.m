function tmatrix = MODtoTOD(julianDate)
%--------------------------------------------------------------------------
% function tmatrix = MODtoTOD(julianDate)
%
% Calculate the transformation matrix from MOD to TOD frame
%--------------------------------------------------------------------------
%
% Inputs:
%     julianDate - Julian date
%     
% Outputs:
%     tmatrix - transformation matrix from MOD (Mean of date) frame to
%               TOD (True of Date) frame
%
% References:     
%     Satellite Orbits Models Methods Applications,
%     Oliver Montenbruck and Eberhard Hill (chapter 5) and Astronomical
%     Alamanac Explanatory Supplement (p 115)
%
%
% Date         Author                   Notes
% ----         ------                   -----
% 06/12/2010   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

TORAD = pi / 180.0;

dpsi = calcNutLon(julianDate) * TORAD;     % nutation in obliquity
eps = calcEclMeanObl(julianDate) * TORAD;  % mean obliquity of the ecliptic
true_eps = calcEclObl(julianDate) * TORAD; % true obliquity of the ecliptic

tmatrix(1,1) =  cos(dpsi);
tmatrix(2,1) =  cos(true_eps) * sin(dpsi);
tmatrix(3,1) =  sin(true_eps) * sin(dpsi);

tmatrix(1,2) = -cos(eps) * sin(dpsi);
tmatrix(2,2) =  cos(eps) * cos(true_eps) * cos(dpsi) + sin(eps) * sin(true_eps);
tmatrix(3,2) =  cos(eps) * sin(true_eps) * cos(dpsi) - sin(eps) * cos(true_eps);

tmatrix(1,3) = -sin(eps) * sin(dpsi);
tmatrix(2,3) =  sin(eps) * cos(true_eps) * cos(dpsi) - cos(eps) * sin(true_eps);
tmatrix(3,3) =  sin(eps) * sin(true_eps) * cos(dpsi) + cos(eps) * cos(true_eps);



 

