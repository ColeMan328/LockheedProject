function [deltaVTotal,TOF] = deltaV4(Velocity,Radius1,Radius2,i,semiMajorAxis)
    deltaVs = 4*Velocity*sin(i/2);

    radiusApogee = Radius2;
    semiMajorTransf = (semiMajorAxis + radiusApogee)/2;
    e = (Radius2/semiMajorTransf)-1;
    initialTureAnomaly = acos((semiMajorAxis*(1-e^2)/(semiMajorAxis)-1));
    Et = -mu/(2*semiMajorTransf);
    vt1 = sqrt(2*(mu/norm(Radius1)-Et));
    phi = atan(e*sin(initialTureAnomaly/(1+e*cos(initialTureAnomaly)));
    deltaVn = 2*sqrt(vt1^2+(norm(Velocity))^2-2*vt1*norm(Velocity)*cos(phi));
    deltaVTotal = deltaVn + deltaVs;

    TOF = maneuverTime(initialTrueAnomaly,semiMajor,e);

end
