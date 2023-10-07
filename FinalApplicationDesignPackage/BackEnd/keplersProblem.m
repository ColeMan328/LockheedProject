% This function finds how long the time of flight is and the final true anomaly. It then calls to the updatedVectors function to find the final position of the satellite
function keplersProblem(L, longitude, H, GST, E, TA, mu, a, TOF, i, RAAN, AOP)
fID = fopen('output.txt','a');
    % Calculating Initial Eccentric Anomaly
    Ei = acos((E + cosd(TA))/(1 + E *cosd(TA))); %radians

    % Half Plane Check
    if TA > 180 % TA is in degrees
        Ei = (2*pi) - Ei;
    end

    % Calculating Inital Mean Anomaly 
    Mi = Ei - (E*sin(Ei));

    % Calculating Mean Motion
    n = sqrt(mu/a^3);

    % Calculating Final Mean Anomaly
    Mf = Mi + (n*TOF);

    % Determining how many times it has orbit
    k = floor(Mf/(2*pi));

    % Recalculating Final Mean Anomaly
    Mf = Mf - k*(2*pi);

    % Half Plane Check
    if Mf > 2*pi
        Mf = Mf - k*2*pi;
    end

    % Newton's Iteration Method to determine Ef
    M0 = Mf - E*sin(Mf);
    E0 = Mf + (Mf- M0)/(1-E*cos(Mf));
    Mn = E0 - E*sin(E0);
    En = E0 + (Mf - Mn)/(1-E*cos(E0));
    epsilon = En-E0;

    while abs(epsilon) >= 0.00001
        Mn = En - E*sin(En);
        En_iteration = En + (Mf - Mn)/ (1 - E*cos(En));

        epsilon = En_iteration - En;
        En = En_iteration;
    end

    Ef = En;

    % Calculating Final True Anomaly
    FTA = acos((cos(Ef) - E)/ (1 - E*cos(Ef)))*(180/pi);
    if Ef > pi
        FTA = 360 - FTA;
    end

    fprintf(fID,"After a Time of Flight of %.4f hours \n", (TOF/(60*60)));
    fprintf(fID,"Final True Anomaly = %.4f \n\n", FTA);
    updatedVectors(L, longitude, H, GST, TOF, a, E, FTA, i, RAAN, AOP, mu)

end