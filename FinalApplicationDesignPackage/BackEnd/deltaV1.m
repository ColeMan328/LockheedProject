% This function determines the delta V and TOF of a user specified Hohmann
% transfer maneuver
function [deltaVTotal,TOF] = deltaV1(Velocity, Radius1, Radius2,i)
    deltaVs = 4*Velocity*sind(i/2); % calculating the delta V required for plane change

    at = (Radius1 + Radius2)/2; % calculating transfer otbit semi-major axis
    Et = -398600.5/(2*at);  % calculating transfer orbit mechanical energy
    Vt1 = sqrt(2*((398600.5/Radius1)+Et));  % calculating velocity of transfer orbit
    dV1 = abs(Vt1-Velocity);    % calculating delta V for first transfer
    V2 = sqrt(398600.5/Radius2);    % calculating velocity of second orbit
    Vt2 = sqrt(2*((398600.5/Radius2)+Et));  % calculating velocity of second transfer orbit
    dV2 = abs(Vt2-V2);  % calculating delta V for second transfer

    Vt3 = sqrt(2*((398600.5/Radius2)+Et));  % calculating velocity of transfer orbit
    dV3 = abs(Vt1-Velocity);    % calculating delta V for third transfer
    Vt4 = sqrt(2*((398600.5/Radius2)+Et));  % calculating velocity of second orbit
    dV4 = abs(Vt2-V2);  % calculating delta V for second transfer
    deltaVTotal = dV1 + dV2 + dV3 + dV4 + deltaVs;  % calculating total delta V
    TOF = 2*pi*sqrt((at^3)/398600.5) *(1/3600); % calculating orbit TOF
end
