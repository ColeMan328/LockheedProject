% This function determines the most effecient Hohmann transfer maneuver for
% the satellite. 
function [minDeltaV,corrRadius,corrTOF,minTOF,corrRadius2,corrDeltaV] = deltaV2(Velocity, Radius1, i)
    deltaVs = 4*Velocity*sind(i/2); % calculating the delta V required for plane change
    for n = 1:(Radius1-8378)
        Radius2(n) = n+8378;    % changing transfer orbit radius
        at(n) = (Radius1 + Radius2(n))/2;   % calculating transfer otbit semi-major axis
        Et(n) = -398600.5/(2*at(n));    % calculating transfer orbit mechanical energy
        Vt1(n) = sqrt(2*((398600.5/Radius1)+Et(n)));    % calculating velocity of transfer orbit
        dV1(n) = abs(Vt1(n)-Velocity);  % calculating delta V for first transfer
        V2(n) = sqrt(398600.5/Radius2(n));  % calculating velocity of second orbit
        Vt2(n) = sqrt(2*((398600.5/Radius2(n))+Et(n))); % calculating velocity of second transfer orbit
        dV2(n) = abs(Vt2(n)-V2(n)); % calculating delta V for second transfer 
        
        Vt3(n) = sqrt(2*((398600.5/Radius2(n))+Et(n))); % calculating velocity of transfer orbit
        dV3(n) = abs(Vt1(n)-Velocity);  % calculating delta V for third transfer
        Vt4(n) = sqrt(2*((398600.5/Radius2(n))+Et(n))); % calculating velocity of second orbit
        dV4(n) = abs(Vt2(n)-V2(n)); % calculating delta V for second transfer
        deltaVTotal(n) = dV1(n) + dV2(n) + dV3(n) + dV4(n) + deltaVs;   % calculating total delta V
        TOF(n) = 2*pi*sqrt((at(n)^3)/398600.5) *(1/3600);   % calculating orbit TOF
    end
    minDeltaV = min(deltaVTotal);   % finding minimum delta V
    index = find(deltaVTotal == minDeltaV); % indexing corrolating minimum values for delta V
    corrRadius = Radius2(index);    % finding radius corresponding to minimum delta V
    corrTOF = TOF(index);   % finding TOF corresponding to minimum delta V
    
    minTOF = min(TOF);  % finding minimum TOF
    index2 = find(TOF == minTOF);   % indexing corrolating minimum values for TOF
    corrRadius2 = Radius2(index2);  % finding radius corresponding to minimum TOF
    corrDeltaV = deltaVTotal(index2);   % finding delta V corresponding to minimum TOF
end
