% This function just prints out the final range, azimuth and elevation
function finalPrinting(Az_final, El_final, rho_final)
fID = fopen('output.txt','a');
    fprintf(fID,"Satellite is located at: \n")
    fprintf(fID,"Final Range = %.4f km \n", rho_final)
    fprintf(fID,"Final Azimuth = %.4f degrees \n", Az_final)
    fprintf(fID,"Final Elevation = %.4f degrees \n\n", El_final)

end