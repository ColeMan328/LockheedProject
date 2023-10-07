% This function just prints out the final range, azimuth and elevation
function finalPrinting(Az_final, El_final, rho_final)

    fprintf("Satellite is located at: \n")
    fprintf("Final Range = %.4f km \n", rho_final)
    fprintf("Final Azimuth = %.4f degrees \n", Az_final)
    fprintf("Final Elevation = %.4f degrees \n\n", El_final)

end