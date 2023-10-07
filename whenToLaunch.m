function TOF = whenToLaunch(semiMajor,ui)
    mu = 398600.5;
    n = sqrt(mu/(semiMajor^3));
    TOF1 = (0-ui)/n;
    TOF2 = (180-ui)/n;
    if TOF1 < 0
        TOF1 = 360 + TOF1;
    end
    if TOF2 < 0
        TOF2 = 360 + TOF2;
    end
    if TOF1 > TOF2
        TOF = TOF2/15;
    else
        TOF = TOF1/15;
    end
end
