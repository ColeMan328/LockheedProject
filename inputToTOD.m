function tMatrix = inputToTOD(Year, Month, Day, Hour, Min, Sec)
    Date=[Year Month Day Hour Min Sec];
    jDate = jdFromDate(Date);
    tMatrix = J2000toTOD(jDate);
end
