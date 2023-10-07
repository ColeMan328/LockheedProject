function MGL = calcMGL (mjd, r, v)
%--------------------------------------------------------------------------
% function MGL = calcMGL (mjd, r, v)
%
% Calculate the Mean Geographic longitude (MGL) given the modified
% modified julian date and state vector (radius and velocity)
%--------------------------------------------------------------------------
%
% Inputs:
%     r   - radius vector [km] TOD frame
%     v   - velocity vector [km/s] TOD frame
%
% Outputs:
%     MGL - mean geographic longitude [deg East] : RAAN + ARG + MEAN_ANOM - GHA
%
% References: 
%
%
% Date         Author                   Notes
% ----         ------                   -----
% 05/22/2020   Carmen Perez de la Cruz  Initialization
% 12/14/2021   Carmen Perez de la Cruz  Standarize
%--------------------------------------------------------------------------
%

gha = calcGHA(mjd);

% From r and v determine 
Kep_elem = rv_to_kep(r,v);

% MGL = RAAN + ARG + MEAN_ANOM - GHA
MGL = mod(Kep_elem.MeanLon - gha + 360.0,360.0);    % deg











    

