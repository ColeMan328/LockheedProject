% This function determines the optimal orbit maneuvers For a tangental to
% non tangential burn transfer orbit maneuver based on minimal delta V and minimal TOF. 
function [minDeltaV,corrRadius,corrTOF,minTOF,corrRadius2,corrDeltaV] = deltaV3(Velocity,Radius1,i,semiMajorAxis)
    deltaVs = 2*Velocity*sind(i/2); % Determining delta V for planc change
    mu = 398600.5;  % gravitational constant of earth
    for n = 1:1000
        e(n) = n/1000;  % changing ecentricity
        initialTureAnomaly(n) = acos((semiMajorAxis*(1-e(n)^2)/(semiMajorAxis)-1)); % calculating initial true anomoly of transfer orbit
        radiusApogee(n) = semiMajorAxis*(1+e(n));   % calculating radius of perigee of the transfer orbit
        semiMajorTransf(n) = (semiMajorAxis + radiusApogee(n))/2;   % calculating semi-major axis of the transfer function
        Et(n) = -mu/(2*semiMajorTransf(n)); % calculating mechanical energy of the transfer orbit
        vt1(n) = sqrt(2*(mu/norm(Radius1)-Et(n)));  % calculating velocity of transfer orbit
        phi(n) = atan(e(n)*sin(initialTureAnomaly(n)/(1+e(n)*cos(initialTureAnomaly(n))))); % angle between initial orbit tangent and the tangent of the new orbit
        deltaVn(n) = sqrt((vt1(n))^2+(norm(Velocity))^2-2*vt1(n)*norm(Velocity)*cos(phi(n)));   % determining delta V of maneuver
        deltaVTotal(n) = deltaVn(n) + deltaVs;  % determining total delta V of maneuver
        TOF(n) = maneuverTime(initialTureAnomaly(n),radiusApogee(n),e(n));  % determining TOF of maneuver
    end
    minDeltaV = min(deltaVTotal);   % determining minimum delta V posible
    index1 = find(deltaVTotal == minDeltaV);    % indexing minimum delta V value corelation to matrix
    corrRadius = radiusApogee(index1);  % retrieving radius corresponding to minimum delta V
    corrTOF = TOF(index1);  % retrieving TOF corresponding to minimum delta V
    minTOF = min(TOF);  % determining minimum TOF posible
    index2 = find(TOF == minTOF);   % indexing minimum TOF value corelation to matrix
    corrRadius2 = semiMajorTransf(index2);  % retrieving radius corresponding to minimum TOF
    corrDeltaV = deltaVTotal(index2);   % retrieving delta V corresponding to minimum TOF


end
