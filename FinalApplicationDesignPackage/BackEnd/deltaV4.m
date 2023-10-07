% This function determines the delta V and TOF For a tangental to
% non tangential burn transfer orbit maneuver
function [deltaVTotal,TOF] = deltaV4(Velocity,Radius1,Radius2,i,semiMajorAxis)
    deltaVs = 4*Velocity*sind(i/2); % Determining delta V for planc change

    radiusApogee = Radius2; 
    semiMajorTransf = (semiMajorAxis + radiusApogee)/2; % calculating semi-major axis of the transfer function
    e = (Radius2/semiMajorTransf)-1;    % calculating ecentricity
    initialTureAnomaly = acos((semiMajorAxis*(1-e^2)/(semiMajorAxis)-1));   % calculating initial true anomoly of transfer orbit
    Et = -mu/(2*semiMajorTransf);   % calculating mechanical energy of the transfer orbit
    vt1 = sqrt(2*(mu/norm(Radius1)-Et));    % calculating velocity of transfer orbit
    phi = atan(e*sin(initialTureAnomaly/(1+e*cos(initialTureAnomaly))));    % angle between initial orbit tangent and the tangent of the new orbit
    deltaVn = 2*sqrt(vt1^2+(norm(Velocity))^2-2*vt1*norm(Velocity)*cos(phi));   % determining delta V of maneuver
    deltaVTotal = deltaVn + deltaVs;    % determining total delta V of maneuver

    TOF = maneuverTime(initialTrueAnomaly,semiMajor,e); % determining TOF of maneuver

end
