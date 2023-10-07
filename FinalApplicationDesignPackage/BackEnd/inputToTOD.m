% This function calculates the transformation matrix to convert the given R
% and V vectors to the apropriate time variable
function tMatrix = inputToTOD(Year, Month, Day, Hour, Min, Sec)
    Date=[Year Month Day Hour Min Sec]; % setting vairable for date given in J2000 file
    jDate = jdFromDate(Date);   % finding julian date
    tMatrix = J2000toTOD(jDate);    % transformation matrix for J2000 to true of date
end
