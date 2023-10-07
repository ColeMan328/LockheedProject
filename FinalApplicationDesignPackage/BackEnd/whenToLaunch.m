% This function finds the time of flight for the tangential to
% non-tangential orbit maneuver
function TOF = whenToLaunch(semiMajor,ui)
    mu = 398600.5;  % gravitational constant of earth 
    ui = ui*(pi/180);   % initial true longitude
    n = sqrt(mu/(semiMajor^3)); % mean motion of the satellite
    P = 2*pi*sqrt((semiMajor^3)/mu);    % Period of the satellite
    TOF1 = (0-ui)/n;    % TOF before ascending node
    TOF2 = (pi-ui)/n;   % TOF before descending node
    if TOF1 < 0
        TOF1 = P + TOF1;    % fixing TOF if negative
    end
    if TOF2 < 0
        TOF2 = P + TOF2;    % fixing TOF if negative
    end
    if TOF1 > TOF2
        TOF = TOF2/3600;    % converting TOF to hours
    else
        TOF = TOF1/3600;    % converting TOF to hours
    end
end
