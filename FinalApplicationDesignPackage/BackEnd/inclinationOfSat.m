% This function determines if the inclination is a special case (i.e.
% equatorial)
function inclinationOfSat(i)
    if (i < 0.001 || i == 180)
        fprintf("<strong>%s </strong>\n","Inclination is Equatorial");
    end
end