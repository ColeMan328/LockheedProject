function [minDeltaV,corrRadius,corrTOF,minTOF,corrRadius2,corrDeltaV] = deltaV3(Velocity,Radius1,i,semiMajorAxis)
    deltaVs = 2*Velocity*sin(i/2);
    mu = 398600.5;
    for n = 1:1000
        e(n) = n/1000;
        initialTureAnomaly(n) = acos((semiMajorAxis*(1-e(n)^2)/(semiMajorAxis)-1));
        radiusApogee(n) = semiMajorAxis*(1+e(n));
        semiMajorTransf(n) = (semiMajorAxis + radiusApogee(n))/2;
        Et(n) = -mu/(2*semiMajorTransf(n));
        vt1(n) = sqrt(2*(mu/norm(Radius1)-Et(n)));
        phi(n) = atan(e(n)*sin(initialTureAnomaly(n)/(1+e(n)*cos(initialTureAnomaly(n)))));
        deltaVn(n) = sqrt((vt1(n))^2+(norm(Velocity))^2-2*vt1(n)*norm(Velocity)*cos(phi(n)));
        deltaVTotal(n) = deltaVn(n) + deltaVs;
        TOF(n) = maneuverTime(initialTureAnomaly(n),radiusApogee(n),e(n));
    end
    minDeltaV = min(deltaVTotal);
    index1 = find(deltaVTotal == minDeltaV);
    corrRadius = radiusApogee(index1);
    corrTOF = TOF(index1);
    minTOF = min(TOF);
    index2 = find(TOF == minTOF);
    corrRadius2 = semiMajorTransf(index2);
    corrDeltaV = deltaVTotal(index2);


end
