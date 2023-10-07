% This function determines the R and V vectors if they were not previously
% given. The vectors determined in this function will be passed to the
% COE's calculations function.
function [R_vector,V_vector] = findingVectors (GST, L, longitude, H, rho, El, Az, LST, rho_dot, El_dot, Az_dot, TOF)
    fID = fopen('output.txt','w');
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
    
    fprintf(fID,"R = " + R_vector(1,1) + " I + " + R_vector(1,2) + " J + " + R_vector(1,3) + " K \n");

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
    fprintf(fID,"V = " + V_vector(1,1) + " I + " + V_vector(1,2) + " J + " + V_vector(1,3) + " K \n\n");

    orbitCalculations(GST, L, longitude, H, R_vector, V_vector, TOF)
    fclose(fID);

end