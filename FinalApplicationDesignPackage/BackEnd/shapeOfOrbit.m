% This function determines the shape of the orbit based on the eccetricity of the orbit 
% which was calculated in the COE calculations function.
function shapeOfOrbit(eccentricity)
fID = fopen('output.txt','a');
    if (eccentricity < .001) 
        fprintf(fID,"Orbit is Circular\n");

    elseif (0 < eccentricity && eccentricity < 1)
        if (abs(1 - eccentricity) < 0.001)
            fprintf(fID,"Orbit is Parabolic\n");
        else 
            fprintf(fID,"Orbit is Elliptical\n");
        end
    elseif (eccentricity > 1)
        fprintf(fID,"Orbit is Hyperbolic\n");
    end
end