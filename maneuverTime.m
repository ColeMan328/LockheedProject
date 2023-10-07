function TOF = maneuverTime(initialTrueAnomaly,semiMajor,e)
    mu = 398600.5;
    n = sqrt(mu/(semiMajor^3));
    finalTrueAnomaly = (180-initialTrueAnomaly)+180;
    Ei = acos((e+cos(initialTrueAnomaly))/(1+(e*cos(initialTrueAnomaly))));
    Mi = Ei - asin(Ei);
    Ef = acos((e+cos(finalTrueAnomaly))/(1+(e*cos(finalTrueAnomaly))));
    Mf = Ef - asin(Ef);
    TOF = (Mf-Mi)/n;
    
end
