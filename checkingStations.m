% This function checks to see if the tracking station is within range of
% the satellite
function checkingStations(El_final)
fID = fopen('output.txt','a');
     % Checking if the longitude of the satellite is within 20 degress of
     % tracking station
    if El_final >= 20 && El_final <= 90
        fprintf(fID,"Satellite is within range\n")
    end

end