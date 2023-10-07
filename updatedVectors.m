% This functions finds the updated position and velocity vectors for the
% satellite.
function updatedVectors(L, longitude, H, GST, TOF, a, E, FTA, i, RAAN, AOP, mu)
fID = fopen('output.txt','a');
    R = (a*(1-E^2))/(1 + E*cosd(FTA));
    R_pqw = [R*cosd(FTA); R*sind(FTA); 0];
    V_pqw = sqrt(mu/(a*(1-E^2)))*[-sind(FTA); E+cosd(FTA); 0];

    %Rotation for inclination
    ROT_i = [1 0 0;
              0 cosd(-i) -sind(-i);
              0 sind(-i) cosd(-i)];

    %Rotation for raan
    ROT_RAAN = [cosd(-RAAN) -sind(-RAAN) 0;
                sind(-RAAN) cosd(-RAAN) 0;
                0 0 1];

    %Rotation for omega
    ROT_AOP = [cosd(-AOP) -sind(-AOP) 0;
                sind(-AOP) cosd(-AOP) 0;
                0 0 1];

    T = ROT_AOP*ROT_i*ROT_RAAN;

    %Updating R and V vectors
    R_Vector_2 = T\R_pqw;
    V_Vector_2 = T\V_pqw;

    fprintf(fID,"Updated R and V vectors: \n" )
    fprintf(fID,"\n")
    fprintf(fID,"R = " + R_Vector_2(1) + " I + " + R_Vector_2(2) + " J + " + R_Vector_2(3) + " K \n");
    fprintf(fID,"V = " + V_Vector_2(1) + " I + " + V_Vector_2(2) + " J + " + V_Vector_2(3) + " K \n\n");
    finalTracking(L, longitude, H, GST, TOF, R_Vector_2)
    findingSuitableTrackingStations(GST, TOF, R_Vector_2)

end