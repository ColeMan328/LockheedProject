function data = COE(R, V)
    % magnitudes/constants
    Rmag = norm(R); % km
    Vmag = norm(V); % km/s
    mu = 398600.5; %km^3/s^2
    
    %% Semi-Major Axis - never undefined
    epsilon = (Vmag^2)/2 - mu/Rmag;
    a = -mu/(2*epsilon); % semi-major axis (km)
    
    %% Eccentricity - never undefined
    h = cross(R,V); % specific angular momentum (km^2/s)
    e_vector = (cross(V,h)/mu) - R/Rmag;
    e = norm(e_vector); % eccentricity
    
    %% Inclination - never undefined
    i = acosd(h(1,3)/norm(h)); % inclination (degrees)
    
    %% Right Ascension of Ascending Node (RAAN)
    khat = [0 0 1];
    n_vector = cross(khat, h);
    Omega = acosd(n_vector(1,1)/norm(n_vector));
    % if equitorial
    if (abs(i) < .001) || (abs(i - 180) < .001)
        % undefined
        Omega = NaN;
    elseif n_vector(1,2) < 0 && (Omega < 180)
        % must be b/t 180 and 360
        Omega = 360 - Omega;
    end
    
    %% Argument of Perigee
    omega = acosd(dot(n_vector, e_vector)/(norm(n_vector)*e));
    % if equitorial or circular orbit
    if (abs(i) < .001) || (abs(i - 180) < .001) || (e < .001)
        omega = NaN;
    elseif e_vector(1,3) < 0 && omega < 180
        omega = 360 - omega;
    end
    
    %% True Anomaly
    nu = acosd(dot(e_vector, R)/(e*Rmag));
    % if circular orbit
    if e < .001
        nu = NaN;
    elseif dot(R, V) < 0
        nu = 360 - nu;
    end
    
    %% Alternate Orbital Elements
    % Argument of Latitude
    u = NaN;
    if (e < 0.001) && (i > 0.001)
        u = acosd(dot(n_vector, R)/(norm(n_vector)*Rmag));
        if R(1,3) < 0
            u = 360 - u;
        end
    end
    
    % Longitude of Perigee
    Pi = NaN;
    if (i < 0.001) && (e > 0.001)
        Pi = acosd(e_vector(1,1)/e);
        if e_vector(1,2) < 0
            Pi = 360 - Pi;
        end
    end
    
    % True Longitude
    l = NaN;
    if (i < 0.001) && (e < 0.001)
        l = acosd(R(1,1)/Rmag);
        if (R(1,2) < 0)
            l = 360 - l;
        end
    end
    
    %% Orbit Shape/Type
    Rp = abs(a)*(1-e);
    
    if e < 0.001
        Orbit = 'Circular';
    elseif (abs(e - 1) < 0.001)
        Orbit = 'Parabolic';
    elseif (e > 0 && e < 1)
        Orbit = 'Elliptical';
    elseif (e > 1)
        Orbit = 'Hyperbolic';
    end
    if (Rp < 6378)
        Orbit = append(Orbit, ' and Ballistic Missile');
    end
    if (abs(i) < .001) || (abs(i - 180) < .001)
        Orbit = append(Orbit, ' and Equitorial');
    end

    a = sprintf('%.6f', a);
    e = sprintf('%.6f', e);
    i = sprintf('%.6f', i);
    omega = sprintf('%.6f', omega);
    Omega = sprintf('%.6f', Omega);
    nu = sprintf('%.6f', nu);
    u = sprintf('%.6f', u);
    Pi = sprintf('%.6f', Pi);
    l = sprintf('%.6f', l);
    data = {a, e, i, omega, Omega, nu, u, Pi, l, Orbit};
end