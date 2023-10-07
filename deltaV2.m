function [minDeltaV,corrRadius,corrTOF,minTOF,corrRadius2,corrDeltaV] = deltaV2(Velocity, Radius1, i)
    deltaVs = 4*Velocity*sin(i/2);
    for n = 1:(Radius1-8378)
        Radius2(n) = n+8378;
        at(n) = (Radius1 + Radius2(n))/2;
        Et(n) = -398600.5/(2*at(n));
        Vt1(n) = sqrt(2*((398600.5/Radius1)+Et(n)));
        dV1(n) = abs(Vt1(n)-Velocity);
        V2(n) = sqrt(398600.5/Radius2(n));
        Vt2(n) = sqrt(2*((398600.5/Radius2(n))+Et(n)));
        dV2(n) = abs(Vt2(n)-V2(n));
        
        Vt3(n) = sqrt(2*((398600.5/Radius2(n))+Et(n)));
        dV3(n) = abs(Vt1(n)-Velocity);
        Vt4(n) = sqrt(2*((398600.5/Radius2(n))+Et(n)));
        dV4(n) = abs(Vt2(n)-V2(n));
        deltaVTotal(n) = dV1(n) + dV2(n) + dV3(n) + dV4(n) + deltaVs;
        TOF(n) = 2*pi*sqrt((at(n)^3)/398600.5) *(1/3600);
    end
    minDeltaV = min(deltaVTotal);
    index = find(deltaVTotal == minDeltaV);
    corrRadius = Radius2(index);
    corrTOF = TOF(index);
    
    minTOF = min(TOF);
    index2 = find(TOF == minTOF);
    corrRadius2 = Radius2(index2);
    corrDeltaV = deltaVTotal(index2);

    %fprintf(' deltaVTotal = %.4f \n TOF = %.2f \n', deltaVTotal,TOF)
end
