% Delta And TOF Calculations Part 1 (Homann Transfer)
function deltaV(Velocity, Radius1, Radius2)

    at = (Radius1 + Radius2)/2;
    Et = -398600.5/(2*at);
    Vt1 = sqrt(2*((398600.5/Radius1)+Et));
    dV1 = abs(Vt1-Velocity);
    V2 = sqrt(398600.5/Radius2);
    Vt2 = sqrt(2*((398600.5/Radius2)+Et));
    dV2 = abs(Vt2-V2);
   
    Vt3 = sqrt(2*((398600.5/Radius2)+Et));
    dV3 = abs(Vt1-Velocity);
    Vt4 = sqrt(2*((398600.5/Radius2)+Et));
    dV4 = abs(Vt2-V2);
    deltaVTotal = dV1 + dV2 + dV3 + dV4;
    TOF = 2*pi*sqrt((at^3)/398600.5) *(1/3600);
    fprintf(' deltaVTotal = %.4f \n TOF = %.2f \n', deltaVTotal,TOF)

end