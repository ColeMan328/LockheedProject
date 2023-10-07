function tmatrix = J2000toTOD(julianDate)
%--------------------------------------------------------------------------
% function tmatrix = J2000toTOD(julianDate)
%
% Calculate the transformation matrix from J2000 to TOD frame
%--------------------------------------------------------------------------
%
% Inputs:
%     julianDate - Julian date
%     
% Outputs:
%     tmatrix - transformation matrix from J2000 to TOD frame
%
% References:     
%     Satellite Orbits Models Methods Applications,
%     Oliver Montenbruck and Eberhard Hill (chapter 5)
%
% Date         Author                   Notes
% ----         ------                   -----
% 06/12/2010   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

tmatrix = MODtoTOD(julianDate) * J2000toMOD(julianDate);



 

