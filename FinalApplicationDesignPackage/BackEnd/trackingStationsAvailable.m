% This function just gives a list of available tracking stations in the
% system. It provides the names of these stations along with their
% longitude, latitude, and elevation.
function trackingStationsAvailable(nameStation, latitude, longitude, elevation)
    fID = fopen('output.txt','a');
    latitude = latitude*(180/pi); % converting to degrees
    longitude = longitude*(180/pi); % converting to degrees
    fprintf(fID,"********************************************************** \n");
    fprintf(fID,"Tracking site: %s \n", nameStation)
    fprintf(fID,"Latitude for tracking station: %.4f degrees \n", latitude)
    fprintf(fID,"Longitude for tracking station: %.4f degrees east \n", longitude)
    fprintf(fID,"Altitude above sea level: %.4f km \n\n", elevation)

end