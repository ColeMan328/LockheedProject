function GHA = calcGHA (mjd)
%--------------------------------------------------------------------------
% function GHA = calcGHA (mjd)
%
% Calculate *mean* greenwich hour angle of Aries given an epoch in
% modified Julian date format
%--------------------------------------------------------------------------
%
% Inputs:
%     mjd - modified Julian date
%
% Outputs:
%     GHA - *mean* greenwich hour angle of Aries [deg] [0-360]
%
% References: 
%     Astronomical almanac
%
% Date         Author                   Notes
% ----         ------                   -----
% 03/30/2005   Carmen Perez de la Cruz  Initialization
% 12/14/2021   Carmen Perez de la Cruz  Standarize
%--------------------------------------------------------------------------


JAN2000  = 51544.5;             % Mjd of 1st jan 2000 at 12h

fracday = mjd - floor(mjd);
t0 = (floor(mjd) - JAN2000) / 36525.0; 
t  = (mjd - JAN2000) / 36525.0;

GHA = (24110.54841 + 8640184.812866 * t0 + 0.093104 * t^2 - ...
       6.2E-06 * t^3) / 86400.0 + 1.002737909350795 * fracday;   % days

GHA = GHA * 360.0;
GHA = mod(GHA, 360.0);