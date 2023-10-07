function jDate = jdFromDate(Date)
%--------------------------------------------------------------------------
% function mjdDate = mjdFromDate(Date)
%
% Convert date vector to modified Julian Date
%--------------------------------------------------------------------------
% 
% Inputs:
%     Date  - Date in vector format [yr mon day hr min sec] or [yr mon day]  
% 
% Outputs:
%     jDate - Julian Date
%
% Date         Author                   Notes
% ----         ------                   -----
% 03/30/2005   Carmen Perez de la Cruz  Initialization
% 12/22/2021   Carmen Perez de la Cruz  Standardize function header
%--------------------------------------------------------------------------
%

% J2000 epoch: January 1, 2000, at 12:00 TT
J2000 = datenum(2000,1,1,12,0,0);
jdJ2000 = 2451545.0;

if length(Date) == 6
   DateTime = datenum(Date(1),Date(2),Date(3),Date(4),Date(5),Date(6));
elseif length(Date) == 3
   DateTime = datenum(Date(1),Date(2),Date(3),0,0,0);
else
   error(' * ERROR * Size of input date vector should be 3 or 6 \n');
end

jDate = DateTime - J2000 + jdJ2000;