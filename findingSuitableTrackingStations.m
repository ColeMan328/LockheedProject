% This function has all the tracking station information and calls to
% trackingStationsAvailable and finalTracking functions.
function findingSuitableTrackingStations(GST, TOF, R_Vector_2)
fID = fopen('output.txt','a');
    % Ground Tracking Stations
    fprintf(fID,"Ground Stations Currently in Use: \n")
    fprintf(fID,"************************************************************** \n");
   
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
    fprintf(fID,"************************************************************** \n");
end

