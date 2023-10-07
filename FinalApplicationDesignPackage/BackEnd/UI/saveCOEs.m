function save(filePath, COEs)
    fid = fopen(filePath, 'w');
    fprintf(fid, 'Semi-major Axis = %s km\n', COEs{1});
    fprintf(fid, 'eccentricity = %s\n', COEs{2});
    fprintf(fid, 'inclination = %s degrees\n', COEs{3});
    fprintf(fid, 'Argument of Perigee = %s\n', COEs{4});
    fprintf(fid, 'RAAN = %s\n', COEs{5});
    fprintf(fid, 'True Anomaly = %s\n', COEs{6});
    fprintf(fid, 'Argument of Latitude = %s\n', COEs{7});
    fprintf(fid, 'Longitude of Perigee = %s\n', COEs{8});
    fprintf(fid, 'True Longitude = %s\n', COEs{9});
    fprintf(fid, 'Orbit = %s\n', COEs{10});

end