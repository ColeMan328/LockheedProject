% This function finds the amount of time the transfer orbit maneuver takes
% to complete. THis function ONLY WORKS WITH CIRCULAR ORBIT ASSUMPTION
function TOF = maneuverTime(initialTrueAnomaly,semiMajor,e)
    mu = 398600.5;  % gravitational constant of earth
    P = 2*pi*sqrt((semiMajor^3)/mu);    % period of the satellite's orbit
    n = sqrt(mu/(semiMajor^3)); % mean motion of the satellite
    finalTrueAnomaly = (180-initialTrueAnomaly)+180;    % final true anomaly of the satellite in the transfer orbit, This value can be changed so the satellite can move to a specific position
    Ei = acos((e+cos(initialTrueAnomaly))/(1+(e*cos(initialTrueAnomaly)))); % initial eccentric anomaly
    Mi = Ei - asin(Ei); % initial mean anomaly
    Ef = acos((e+cos(finalTrueAnomaly))/(1+(e*cos(finalTrueAnomaly)))); % initial eccentric anomaly
    Mf = Ef - asin(Ef); % final mean anomaly
    TOF = (Mf-Mi)/n;    % time of flight for orbit maneuver in seconds
    if TOF < 0
        TOF = P + TOF;  % Checking for negative time and corrosponding TOF
    end
    TOF = TOF/3600; % time of flight for orbit maneuver in hours
end
