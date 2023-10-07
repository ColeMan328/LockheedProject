% This function updates the final tracking station and calls to the
% finalPrinting function.
function finalTracking(L, longitude, H, GST, TOF, R_Vector_2)
fID = fopen('output.txt','a');
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
    fprintf(fID,"********************************************************** \n");
end