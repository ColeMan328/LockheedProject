% Taylor Tallerday

clear; clc;

%Case 1:
fprintf("\n <strong>%s </strong>\n","CASE 1:");
fprintf("************************************************************** \n");
    %Site location:
        GST1 = (17+(40/60))*15*(pi/180); %Greenwich Sidereal Time converted to military and in radians
        longitude1 = -64.7*(pi/180); %radians (degrees were given in west therfore negative)
        L1 = 76.53*(pi/180); % (latitude in radians)
        H1 = 1207/1000; %km
        LST1 = GST1 + longitude1; %radians
    %Satellite location:
        rho1 = 3000; %km
        Az1 = 7.5*(pi/180); %radians/s
        El1 = 85*(pi/180); %radians
        rho_dot1 = 6; %km/s
        Az_dot1 = 1*(pi/180); %rad/s
        El_dot1 = 0.01*(pi/180); %rad/s
        %Time of Flight 2 weeks
        TOF = 2*7*24*3600; % (Time of Flight in seconds)
    % Calling to functions
    displayForSiteAndSatellite (GST1, longitude1, H1, rho1, El1, Az1, L1, LST1, rho_dot1, El_dot1, Az_dot1);
    findingVectors(GST1, L1, longitude1, H1, rho1, El1, Az1, LST1, rho_dot1, El_dot1, Az_dot1, TOF);

% Case 2:
fprintf("\n <strong>%s </strong>\n","CASE 2:");
fprintf("************************************************************** \n");
    %Site location:
        GST2 = 22*15*(pi/180); %Greenwich Sidereal Time converted to military and in radians
        longitude2 = (-104.54)*(pi/180); %radians (degrees were given in west therefore negative)
        L2 = (38.8)*(pi/180); % (latitude in radians)
        H2 = 1915/1000; %km
        LST2 = GST2 + longitude2; %radians
    %Satellite location:
        rho2 = 2121.4180; %km
        Az2 = 350*(pi/180); %radians/s
        El2 = 35.3507*(pi/180); %radians
        rho_dot2 = -3.32040; %km/s
        Az_dot2 = -0.07653*(pi/180); %rad/s
        El_dot2 = 0.20367*(pi/180); %rad/s
        %Time of Flight 4 days
        TOF2 = 4*24*3600; % (Time of Flight in seconds)
    % Calling to functions
    displayForSiteAndSatellite (GST2, longitude2, H2, rho2, El2, Az2, L2, LST2, rho_dot2, El_dot2, Az_dot2);
    findingVectors(GST2, L2, longitude2, H2, rho2, El2, Az2, LST2, rho_dot2, El_dot2, Az_dot2, TOF2);
   
% COE Calculations functions. R and V vectors are passed into this function and altitude,
% specific angular momentum, inclination, RAAN, AOP, true anomaly, and
% alternative COE's are all calculated in this function. Plane checks are accounted for
% in these calcuations as well. Calls to other functions that determine the shape of the
% orbit and the special inclination cases.
function orbitCalculations(GST, L, longitude, H, R_vector, V_vector, TOF)
disp("COE's for the spacecraft are:");
fprintf("************************************************************** \n");
%Establishing Local Constants
mu = 398600.5; % km^3/s^2
k_vector = [0, 0, 1];
i_vector = [1, 0, 0];

        % Finding magnitudes for R and V
        R = norm(R_vector);
        V = norm(V_vector);

        % Calculating Altitude
        epsilon = ((V^2)/2) - (mu/R); % Specific Mechanical Energy
        fprintf("Specific Mechanical Energy (epsilon): %.4f\n", epsilon);
       
        a = -mu/(2*epsilon);
        fprintf("Semi-Major Axis: %.4f\n", a);

        %Calculating eccentricity using specific angular momentum
        h_vector = cross(R_vector, V_vector);
        E_vector = ((cross(V_vector, h_vector))/mu) - (R_vector/R);
        E = norm(E_vector);
        fprintf("Eccentricity: %.4f\n", E);

        % Calculating inclination
        h = norm(h_vector);
        i = acosd(h_vector(1,3)/h);
        fprintf("Inclination: %.4f\n", i);

        %Determining if Equatorial
        inclinationOfSat(i);
         n_vector = cross(k_vector, h_vector); % Node Vector
         n = norm(n_vector);

        %Right Ascension of Ascending Node and Argument of Perigee
        if i == 0 || i == 180
            %When it doesn't exist
            fprintf("Right Ascension of Ascending Node: NaN \n");
            fprintf("Argument of Perigee: NaN \n");
        else

            %Calculating Right Ascension of Ascending Node
            RAAN = acosd(n_vector(1,1)/ n);
            if n_vector(1,2) < 0
                RAAN = 360 - RAAN;
            end
            fprintf("Right Ascension of Ascending Node: %.4f\n", RAAN);
            if -0.001 < E && E < 0.001
                %When it doesn't exist
                fprintf("Argument of Perigee: NaN \n");
            else

                %Calculating Arguement of Perigee
                AOP = acosd(dot(n_vector,E_vector)/(n*E));
                if E_vector(1,3) < 0
                    AOP = 360 - AOP;
                end
                fprintf("Argument of Perigee: %.4f\n", AOP);
            end
        end

        %True Anomaly
        if -0.001 < E && E < 0.001
            %When it doesn't exist
            fprintf("True Anomaly: NaN \n");
        else

            %Calculating True Anomaly
            TA = acosd(dot(E_vector,R_vector)/(E*R));
            if dot(R_vector, V_vector) < 0
                TA = 360 - TA;
            end
            fprintf("True Anomaly: %.4f\n", TA);
        end

    %Alternative COE's        
    if i < 0.001
        if -0.001 < E && E < 0.001

            %Calculating True Longitude
            TL = acosd(dot(i_vector, R_vector)/ R);
            if E_vector(1,2) < 0
                TL = 360 - TL;
            end
            fprintf("True Longitude: %.4f\n", TL);
        else

        %Calculating Longitude of Perigee
        LP = acosd(dot(i_vector,E_vector)/E);
        if E_vector(1,2) < 0
                LP = 360 - LP;
        end
        fprintf("Longitude of Perigee: %.4f\n", LP);
        end
    else
        if E < 0.001
            
            %Argument of Latitude
            AL = acosd(dot(n_vector,R_vector)/(n*R));
            if R_vector(1,3) < 0
                AL = 360 - AL;
            end
            fprintf("Argument of Latitude: %.4f\n", AL);
        end
    end

        %Radius of Perigee
        RP = abs(a*(1-E));
        if (RP < 6628)
            fprintf("Ballistic Missile with a Radius of Perigee: %.4f \n", RP);
        end

        disp("**************************************************************");

        %Determining Shape of Orbit
        shapeOfOrbit(E);
        disp("**************************************************************");

        % Calling to Kepler's Problem Function
        keplersProblem(L, longitude, H, GST, E, TA, mu, a, TOF, i, RAAN, AOP)

end

% This function determines the shape of the orbit based on the eccetricity of the orbit
% which was calculated in the COE calculations function.
function shapeOfOrbit(eccentricity)

    if (eccentricity < .001)
        fprintf("<strong>%s </strong>\n","Orbit is Circular");

    elseif (0 < eccentricity && eccentricity < 1)

        if (abs(1 - eccentricity) < 0.001)
            fprintf("<strong>%s </strong>\n","Orbit is Parabolic");
        else
            fprintf("<strong>%s </strong>\n","Orbit is Elliptical");
        end

    elseif (eccentricity > 1)
        fprintf("<strong>%s </strong>\n","Orbit is Hyperbolic");
    end

end

% This function determines if the inclination is a special case (i.e.
% equatorial)
function inclinationOfSat(i)

    if (i < 0.001 || i == 180)
        fprintf("<strong>%s </strong>\n","Inclination is Equatorial");
    end

end

% This function determines the R and V vectors if they were not previously
% given. The vectors determined in this function will be passed to the
% COE's calculations function.
function findingVectors (GST, L, longitude, H, rho, El, Az, LST, rho_dot, El_dot, Az_dot, TOF)

    % Constants
    w_vector = [0; 0; (15*(pi/(180*3600)))]; %Earth's velocity vector in rad/s
    a_earth = 6378.137;
    e_earth = 0.08182;

    % Calculations to determine R vector
    % Finding rho (position) vector
    rho_s = -rho*cos(El)*cos(Az); % S components
    rho_e = rho*cos(El)*sin(Az); % E component
    rho_z = rho*sin(El); % Z component

    rho_SEZ = [rho_s; rho_e; rho_z]; % Position Vector in the SEZ coordinate frame
    % Transform P vector from SEZ to IJK
    rho_IJK = [sin(L)*cos(LST), -sin(LST), cos(L)*cos(LST);
             sin(L)*sin(LST), cos(LST), cos(L)*sin(LST);
             -cos(L),           0,      sin(L)] * rho_SEZ;

    % Determining the position of the site
    % Use long expressions
    x = ((a_earth)/(sqrt(1-((e_earth^2)*sin(L)^2))) + H) * cos(L);
    z = ((a_earth*(1-e_earth^2))/(sqrt(1-((e_earth^2)*sin(L)^2))) + H) * sin(L);
    R_site = [x*cos(LST); x*sin(LST); z];

    % Determining the R vector that will be passed to the calculations
    % function. The vector has to be transposed to allow it to be a vector.
    R_vector = R_site + rho_IJK;
    R_vector = R_vector';
   
    fprintf("R = " + R_vector(1,1) + " I + " + R_vector(1,2) + " J + " + R_vector(1,3) + " K \n");

    % Calculations to determine V vector
    rho_dot_s = -rho_dot*cos(El)*cos(Az) + rho*El_dot*sin(El)*cos(Az) + rho*Az_dot*cos(El)*sin(Az);
    rho_dot_e = rho_dot*cos(El)*sin(Az) - rho*El_dot*sin(El)*sin(Az) + rho*Az_dot*cos(El)*cos(Az);
    rho_dot_z = rho_dot*sin(El) + rho*El_dot*cos(El);

    rho_dot_SEZ = [rho_dot_s; rho_dot_e; rho_dot_z];

    rho_dot_IJK = [sin(L)*cos(LST), -sin(LST), cos(L)*cos(LST);
             sin(L)*sin(LST), cos(LST), cos(L)*sin(LST);
             -cos(L),           0,      sin(L)]* rho_dot_SEZ;

    rho_dot_IJK = rho_dot_IJK';
    w_R_cross = cross(w_vector, R_vector);

    % Determining the V vector that will be passed to the calculations
    % function.The vector has to be transposed to allow it to be a vector.
    V_vector = rho_dot_IJK + w_R_cross;
    fprintf("V = " + V_vector(1,1) + " I + " + V_vector(1,2) + " J + " + V_vector(1,3) + " K \n\n");
   
    % Calling to the function to make calculations for the orbit
    orbitCalculations(GST, L, longitude, H, R_vector, V_vector, TOF);

end

% This function just displays the site and satellite location values that
% were used in the calculations in a nicely formatted way.
function displayForSiteAndSatellite (GST, longitude, H, rho, El, Az, L, LST, rho_dot, El_dot, Az_dot)

% Converting all units to degrees
L = L/(pi/180);
longitude = longitude/(pi/180);
GST = GST/(pi/180);
LST = LST/(pi/180);
Az = Az/(pi/180);
El = El/(pi/180);

    fprintf("Inital Range tracking site location is at: \n")
    fprintf("Latitude = %.4f degrees \n", L)
    fprintf("Longitude = %.4f degrees east \n", longitude)
    fprintf("Altitude above sea level = %.4f km \n", H)
    fprintf("Greenwich Sidereal Time = %.4f degrees \n", GST)
    fprintf("Local Sidereal Time = %.4f degrees \n\n", LST)
    fprintf("Satellite is located at: \n")
    fprintf("Range = %.4f km \n", rho)
    fprintf("Azimuth = %.4f degrees \n", Az)
    fprintf("Elevation = %.4f degrees \n\n", El)
    fprintf("Satellite is moving at: \n")
    fprintf("Range rate of change = %.4f km/sec\n", rho_dot)
    fprintf("Azimuth rate of change = %.4f rad/sec\n", Az_dot)
    fprintf("Elevation rate of change = %.4f rad/sec\n\n\n", El_dot)

end

% This function finds how long the time of flight is and the final true
% anomaly. It then calls to the updatedVectors function to find the final
% position of the satellite
function keplersProblem(L, longitude, H, GST, E, TA, mu, a, TOF, i, RAAN, AOP)

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

    fprintf("After a Time of Flight of %.4f hours \n", (TOF/(60*60)));
    fprintf("Final True Anomaly = %.4f \n\n", FTA);
    % Calling to function to find updated R and V vectors
    updatedVectors(L, longitude, H, GST, TOF, a, E, FTA, i, RAAN, AOP, mu)

end

% This functions finds the updated position and velocity vectors for the
% satellite.
function updatedVectors(L, longitude, H, GST, TOF, a, E, FTA, i, RAAN, AOP, mu)
    R = (a*(1-E^2))/(1 + E*cosd(FTA));
    R_pqw = [R*cosd(FTA); R*sind(FTA); 0];
    V_pqw = sqrt(mu/(a*(1-E^2)))*[-sind(FTA); E+cosd(FTA); 0];

    %Rotation for Inclination
    ROT_i = [1 0 0;
              0 cosd(-i) -sind(-i);
              0 sind(-i) cosd(-i)];

    %Rotation for Raan
    ROT_RAAN = [cosd(-RAAN) -sind(-RAAN) 0;
                sind(-RAAN) cosd(-RAAN) 0;
                0 0 1];

    %Rotation for Argument of Perigee
    ROT_AOP = [cosd(-AOP) -sind(-AOP) 0;
                sind(-AOP) cosd(-AOP) 0;
                0 0 1];
    T = ROT_AOP*ROT_i*ROT_RAAN;

    %Updating R and V vectors
    R_Vector_2 = T\R_pqw;
    V_Vector_2 = T\V_pqw;

    fprintf("Final Range Tracking Site Location at: \n")
    fprintf("R = " + R_Vector_2(1) + " I + " + R_Vector_2(2) + " J + " + R_Vector_2(3) + " K \n");
    fprintf("V = " + V_Vector_2(1) + " I + " + V_Vector_2(2) + " J + " + V_Vector_2(3) + " K \n\n");  
    finalTracking(L, longitude, H, GST, TOF, R_Vector_2)
    findingSuitableTrackingStation(GST, TOF, R_Vector_2)

    %deltaV(sqrt(sum(V_Vector_2.^2)),sqrt(sum(R_Vector_2.^2)),20000)

end

% This function updates the final tracking station and calls to the
% finalPrinting function.
function finalTracking(L, longitude, H, GST, TOF, R_Vector_2)

    % Converting units
    L = L/(pi/180); % rad to degrees
    longitude = longitude/(pi/180); % rad to degrees
    GST = GST/(pi/180); % rad to degrees
    TOF = TOF/3600; % seconds to hours
       
    % Constants
    a_earth = 6378.137;
    e_earth = 0.08182;

% Calculating Final Tracking Station Variables
    GST_final = GST + TOF*15;
    LST_final = GST_final + longitude;

    % Determining the position of the site
    % Use long expressions
    x = ((a_earth)/(sqrt(1-((e_earth^2)*sind(L)^2))) + H) * cosd(L); 
    z = ((a_earth*(1-e_earth^2))/(sqrt(1-((e_earth^2)*sind(L)^2))) + H) * sind(L);
    R_site_new = [x*cosd(LST_final); x*sind(LST_final); z];
    rho_IJK_new = R_Vector_2 - R_site_new;
    ROT = [sind(L)*cosd(LST_final) -sind(LST_final) cosd(L)*cosd(LST_final);
           sind(L)*sind(LST_final) cosd(LST_final) cosd(L)*sind(LST_final);
           -cosd(L)           0      sind(L)];
    rho_dot_SEZ_new = ROT\rho_IJK_new;

    % Calculating the range, elevation, and azimuth
    rho_final = norm(rho_dot_SEZ_new);
    El_final = asind(rho_dot_SEZ_new(3)/rho_final);
    Az_final = acosd(-rho_dot_SEZ_new(1)/(rho_final*cosd(El_final)));

    % Plane Check
    if rho_dot_SEZ_new(2) < 0
        Az_final = 360 - Az_final;
    end

    finalPrinting(Az_final, El_final, rho_final)
    checkingStations(El_final)
end

% This function just prints out the final range, azimuth and elevation
function finalPrinting(Az_final, El_final, rho_final)
    fprintf("Satellite is located at: \n")
    fprintf("Final Range = %.4f km \n", rho_final)
    fprintf("Final Azimuth = %.4f degrees \n", Az_final)
    fprintf("Final Elevation = %.4f degrees \n\n", El_final)
end

% This function just gives a list of available tracking stations in the
% system. It provides the names of these stations along with their
% longitude, latitude, and elevation.
function trackingStationsAvailable(nameStation, latitude, longitude, elevation)

    fprintf("Tracking site: %s \n", nameStation)
    fprintf("Latitude for tracking station: %.4f \n", latitude)
    fprintf("Longitude for tracking station: %.4f \n", longitude)
    fprintf("Altitude above sea level: %.4f \n\n", elevation)

end

% This function has all the tracking station information and calls to
% trackingStationsAvailable and finalTracking functions.
function findingSuitableTrackingStation(GST, TOF, R_Vector_2)

    % Ground Tracking Stations
    fprintf("<strong>%s </strong>\n","Ground Stations Currently in Use:")
    fprintf("************************************************************** \n");
   
    % Colorado (CTS)
    nameStation_CTS = "Colorado";
    latitude_CTS = 38.8199*pi/180; % degrees N
    longitude_CTS = -104.7004*pi/180; % degrees E
    elevation_CTS = 1875/1000; % km
    trackingStationsAvailable(nameStation_CTS, latitude_CTS, longitude_CTS, elevation_CTS)
    finalTracking(latitude_CTS, longitude_CTS, elevation_CTS, GST, TOF, R_Vector_2)
   
    % Diego Garcia (DGS)
    nameStation_DGS = "Diego Garcia";
    latitude_DGS = -7.3192*pi/180; % degrees N
    longitude_DGS = 72.4227*pi/180; % degrees E
    elevation_DGS = 3/1000; % km
    trackingStationsAvailable(nameStation_DGS, latitude_DGS, longitude_DGS, elevation_DGS)
    finalTracking(latitude_DGS, longitude_DGS, elevation_DGS, GST, TOF, R_Vector_2)
   
    % Guam (GTS)
    nameStation_GTS = "Guam";
    latitude_GTS = 13.308610*pi/180; % degrees N
    longitude_GTS = -144.734440*pi/180; % degrees E
    elevation_GTS = 0.082000; % km
    trackingStationsAvailable(nameStation_GTS, latitude_GTS, longitude_GTS, elevation_GTS)
    finalTracking(latitude_GTS, longitude_GTS, elevation_GTS, GST, TOF, R_Vector_2)
   
    % Hawaii (HTS)
    nameStation_HTS = "Hawaii";
    latitude_HTS = 21.5577*pi/180; % degrees N
    longitude_HTS = -158.2487*pi/180; % degrees E
    elevation_HTS = 132/1000; % km
    trackingStationsAvailable(nameStation_HTS, latitude_HTS, longitude_HTS, elevation_HTS)
    finalTracking(latitude_HTS, longitude_HTS, elevation_HTS, GST, TOF, R_Vector_2)

    % New Hampshire (NHS)
    nameStation_NHS = "New Hampshire";
    latitude_NHS = 42.9804*pi/180; % degrees N
    longitude_NHS = -71.6828*pi/180; % degrees E
    elevation_NHS = 335/1000; % km
    trackingStationsAvailable(nameStation_NHS, latitude_NHS, longitude_NHS, elevation_NHS)
    finalTracking(latitude_NHS, longitude_NHS, elevation_NHS, GST, TOF, R_Vector_2)
   
    %Telemetry & Command Station
    nameStation_REF = "Telemetry & Command Station";
    latitude_REF = 51.1319*pi/180; % degrees N
    longitude_REF = -0.8963*pi/180; % degrees E
    elevation_REF = 147/1000; % km
    trackingStationsAvailable(nameStation_REF, latitude_REF, longitude_REF, elevation_REF)
    finalTracking(latitude_REF, longitude_REF, elevation_REF, GST, TOF, R_Vector_2)

    % Thule (TTS)
    nameStation_TTS = "Thule";
    latitude_TTS = 76.5319*pi/180; % degrees N
    longitude_TTS = -68.7032*pi/180; % degrees E
    elevation_TTS = 76/1000; % km
    trackingStationsAvailable(nameStation_TTS, latitude_TTS, longitude_TTS, elevation_TTS)
    finalTracking(latitude_TTS, longitude_TTS, elevation_TTS, GST, TOF, R_Vector_2)

    % Vandenberg (VTS)
    nameStation_VTS = "Vandenberg";
    latitude_VTS = 34.6348*pi/180; % degrees N
    longitude_VTS = -120.6122*pi/180; % degrees E
    elevation_VTS = 153/1000; % km
    trackingStationsAvailable(nameStation_VTS, latitude_VTS, longitude_VTS, elevation_VTS)
    finalTracking(latitude_VTS, longitude_VTS, elevation_VTS, GST, TOF, R_Vector_2)
    fprintf("************************************************************** \n");
end

% This function checks to see if the tracking station is within range of
% the satellite
function checkingStations(El_final)

     % Checking if the longitude of the satellite is within 20 degress of
     % tracking station
    if El_final >= 20 && El_final <= 90
        fprintf("<strong>%s </strong>\n","Satellite is within range")
    end

end

% Delta And TOF Calculations Part 1 (Homann Transfer)
function deltaV(Velocity, Radius1, Radius2)

    at = (Radius1 + Radius2)/2;
    Et = -398600.5/(2*at);
    Vt1 = sqrt(2*((398600.5/Radius1)+Et));
    dV1 = abs(Vt1-Velocity);
    V2 = sqrt(398600.5/Radius2);
    Vt2 = sqrt(2*((398600.5/Radius2)+Et));
    dV2 = abs(Vt2-V2);
   
    Vt3 = sqrt(2*((398600.5/Radius2)+Et));
    dV3 = abs(Vt1-Velocity);
    Vt4 = sqrt(2*((398600.5/Radius2)+Et));
    dV4 = abs(Vt2-V2);
    deltaVTotal = dV1 + dV2 + dV3 + dV4;
    TOF = 2*pi*sqrt((at^3)/398600.5) *(1/3600);
    fprintf(' deltaVTotal = %.4f \n TOF = %.2f \n', deltaVTotal,TOF)

end

% Function J2000toTOD()
%}