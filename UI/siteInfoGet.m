function siteData = getSiteInfo(R, TOF, GST)
    GST = deg2rad((GST*15));
    R = R';
    [site, lat, long, alt, rho, Az, El] = findingSuitableTrackingStation(GST, TOF, R);
    lat = sprintf("%.6f", lat);
    long = sprintf("%.6f", long);
    alt = sprintf("%.6f", alt);
    rho = sprintf("%.6f", rho);
    Az = sprintf("%.6f", Az);
    El = sprintf("%.6f", El);
    siteData = {site, lat, long, alt, rho, Az, El};

    % This function updates the final tracking station and calls to the
    % finalPrinting function.
    function [Az_final, El_final, rho_final] = finalTracking(L, longitude, H, GST, TOF, R_Vector_2)
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
        x = ((a_earth)/(sqrt(1-((e_earth^2)*sind(L)^2))) + H) * cosd(L); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    end

    % This function checks to see if the tracking station is within range of
    % the satellite
    function inRange = checkingStations(Az_final, El_final)
        % Checking if the longitude of the satellite is within 20 degress of
        % tracking station
        if El_final >= 20 && El_final <= 90
            inRange = 1;
        else
            inRange = 0;
        end
    end

    % This function has all the tracking station information and calls to
    % trackingStationsAvailable and finalTracking functions.
    function [site, lat, long, alt, rho, Az, El] = findingSuitableTrackingStation(GST, TOF, R_Vector_2)
        % Ground Tracking Stations
            
        % Colorado (CTS)
        nameStation_CTS = "Colorado";
        latitude_CTS = 38.8199*pi/180; % degrees N
        longitude_CTS = -104.7004*pi/180; % degrees E
        elevation_CTS = 1875/1000; % km
        [az_CTS, el_CTS, rho_CTS] = finalTracking(latitude_CTS, longitude_CTS, elevation_CTS, GST, TOF, R_Vector_2);
        if (checkingStations(az_CTS, el_CTS)) == 1
            site = nameStation_CTS;
            lat = latitude_CTS;
            long = longitude_CTS;
            alt = elevation_CTS;
            rho = rho_CTS;
            Az = az_CTS;
            El = el_CTS;
        end
        % Diego Garcia (DGS)
        nameStation_DGS = "Diego Garcia";
        latitude_DGS = -7.3192*pi/180; % degrees N
        longitude_DGS = 72.4227*pi/180; % degrees E
        elevation_DGS = 3/1000; % km
        [az_DGS, el_DGS, rho_DGS] = finalTracking(latitude_DGS, longitude_DGS, elevation_DGS, GST, TOF, R_Vector_2);
        if (checkingStations(az_DGS, el_DGS)) == 1
            site = nameStation_DGS;
            lat = latitude_DGS;
            long = longitude_DGS;
            alt = elevation_DGS;
            rho = rho_DGS;
            Az = az_DGS;
            El = el_DGS;
        end

        % Guam (GTS)
        nameStation_GTS = "Guam";
        latitude_GTS = 13.308610*pi/180; % degrees N
        longitude_GTS = -144.734440*pi/180; % degrees E
        elevation_GTS = 0.082000; % km
        [az_GTS, el_GTS, rho_GTS] = finalTracking(latitude_GTS, longitude_GTS, elevation_GTS, GST, TOF, R_Vector_2);
        if (checkingStations(az_GTS, el_GTS)) == 1
            site = nameStation_GTS;
            lat = latitude_GTS;
            long = longitude_GTS;
            alt = elevation_GTS;
            rho = rho_GTS;
            Az = az_GTS;
            El = el_GTS;
        end

        % Hawaii (HTS)
        nameStation_HTS = "Hawaii";
        latitude_HTS = 21.5577*pi/180; % degrees N
        longitude_HTS = -158.2487*pi/180; % degrees E
        elevation_HTS = 132/1000; % km
        [az_HTS, el_HTS, rho_HTS] = finalTracking(latitude_HTS, longitude_HTS, elevation_HTS, GST, TOF, R_Vector_2);
        if (checkingStations(az_HTS, el_HTS)) == 1
            site = nameStation_HTS;
            lat = latitude_HTS;
            long = longitude_HTS;
            alt = elevation_HTS;
            rho = rho_HTS;
            Az = az_HTS;
            El = el_HTS;
        end

        % New Hampshire (NHS)
        nameStation_NHS = "New Hampshire";
        latitude_NHS = 42.9804*pi/180; % degrees N
        longitude_NHS = -71.6828*pi/180; % degrees E
        elevation_NHS = 335/1000; % km
        [az_NHS, el_NHS, rho_NHS] = finalTracking(latitude_NHS, longitude_NHS, elevation_NHS, GST, TOF, R_Vector_2);
        if (checkingStations(az_NHS, el_NHS)) == 1
            site = nameStation_NHS;
            lat = latitude_NHS;
            long = longitude_NHS;
            alt = elevation_NHS;
            rho = rho_NHS;
            Az = az_NHS;
            El = el_NHS;
        end
   
        %Telemetry & Command Station
        nameStation_REF = "Telemetry & Command Station";
        latitude_REF = 51.1319*pi/180; % degrees N
        longitude_REF = -0.8963*pi/180; % degrees E
        elevation_REF = 147/1000; % km
        [az_REF, el_REF, rho_REF] = finalTracking(latitude_REF, longitude_REF, elevation_REF, GST, TOF, R_Vector_2);
        if (checkingStations(az_REF, el_REF)) == 1
            site = nameStation_REF;
            lat = latitude_REF;
            long = longitude_REF;
            alt = elevation_REF;
            rho = rho_REF;
            Az = az_REF;
            El = el_REF;
        end

        % Thule (TTS)
        nameStation_TTS = "Thule";
        latitude_TTS = 76.5319*pi/180; % degrees N
        longitude_TTS = -68.7032*pi/180; % degrees E
        elevation_TTS = 76/1000; % km
        [az_TTS, el_TTS, rho_TTS] = finalTracking(latitude_TTS, longitude_TTS, elevation_TTS, GST, TOF, R_Vector_2);
        if (checkingStations(az_TTS, el_TTS)) == 1
            site = nameStation_TTS;
            lat = latitude_TTS;
            long = longitude_TTS;
            alt = elevation_TTS;
            rho = rho_TTS;
            Az = az_TTS;
            El = el_TTS;
        end

        % Vandenberg (VTS)
        nameStation_VTS = "Vandenberg";
        latitude_VTS = 34.6348*pi/180; % degrees N
        longitude_VTS = -120.6122*pi/180; % degrees E
        elevation_VTS = 153/1000; % km
        [az_VTS, el_VTS, rho_VTS] = finalTracking(latitude_VTS, longitude_VTS, elevation_VTS, GST, TOF, R_Vector_2);
        if (checkingStations(az_VTS, el_VTS)) == 1
            site = nameStation_VTS;
            lat = latitude_VTS;
            long = longitude_VTS;
            alt = elevation_VTS;
            rho = rho_VTS;
            Az = az_VTS;
            El = el_VTS;
        end
    end
end