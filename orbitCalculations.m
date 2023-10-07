% COE Calculations functions. R and V vectors are passed into this function and altitude,
% specific angular momentum, inclination, RAAN, AOP, true anomaly, and
% alternative COE's are all calculated in this function. Plane checks are accounted for 
% in these calcuations as well. Calls to other functions that determine the shape of the 
% orbit and the special inclination cases.
function orbitCalculations(GST, L, longitude, H, R_vector, V_vector, TOF)
fID = fopen('output.txt','a');

fprintf(fID,"COE's for the spacecraft are:\n");
fprintf(fID,"************************************************************** \n");

%Establishing Local Constants
mu = 398600.5; % km^3/s^2
k_vector = [0, 0, 1];
i_vector = [1, 0, 0];

        % Finding magnitudes for R and V
        R = norm(R_vector);
        V = norm(V_vector);

        % Calculating Altitude
        epsilon = ((V^2)/2) - (mu/R); % Specific Mechanical Energy
        fprintf(fID,"Specific Mechanical Energy (epsilon): %.4f\n", epsilon);
        
        a = -mu/(2*epsilon);
        fprintf(fID,"Semi-Major Axis: %.4f\n", a);

        %Calculating eccentricity using specific angular momentum
        h_vector = cross(R_vector, V_vector);
        E_vector = ((cross(V_vector, h_vector))/mu) - (R_vector/R);
        E = norm(E_vector);
        fprintf(fID,"Eccentricity: %.4f\n", E);

        % Calculating inclination 
        h = norm(h_vector);
        i = acosd(h_vector(1,3)/h);
        fprintf(fID,"Inclination: %.4f\n", i);

        %Determining if Equatorial
        inclinationOfSat(i);

         n_vector = cross(k_vector, h_vector); % Node Vector
         n = norm(n_vector);

        %Right Ascension of Ascending Node and Argument of Perigee
        if i == 0 || i == 180
            %When it doesn't exist
            fprintf(fID,"Right Ascension of Ascending Node: NaN \n");
            fprintf(fID,"Argument of Perigee: NaN \n");
        else
            %Calculating Right Ascension of Ascending Node
            RAAN = acosd(n_vector(1,1)/ n);
            if n_vector(1,2) < 0
                RAAN = 360 - RAAN;
            end
            fprintf(fID,"Right Ascension of Ascending Node: %.4f\n", RAAN);

            if -0.001 < E && E < 0.001
                %When it doesn't exist
                AOP = [];
                fprintf(fID,"Argument of Perigee: NaN \n");
            else
                %Calculating Arguement of Perigee
                AOP = acosd(dot(n_vector,E_vector)/(n*E));
                if E_vector(1,3) < 0
                    AOP = 360 - AOP;
                end
                fprintf(fID,"Argument of Perigee: %.4f\n", AOP);
            end
        end

        %True Anomaly
        if -0.001 < E && E < 0.001
            %When it doesn't exist
            TA = [];
            fprintf(fID,"True Anomaly: NaN \n");
        else 
            %Calculating True Anomaly
            TA = acosd(dot(E_vector,R_vector)/(E*R));
            if dot(R_vector, V_vector) < 0
                TA = 360 - TA;
            end
            fprintf(fID,"True Anomaly: %.4f\n", TA);
        end

    %Alternative COE's        
    if i < 0.001
        if -0.001 < E && E < 0.001
            %Calculating True Longitude
            TL = acosd(dot(i_vector, R_vector)/ R);
            if E_vector(1,2) < 0
                TL = 360 - TL;
            end
            fprintf(fID,"True Longitude: %.4f\n", TL);
        else
        %Calculating Longitude of Perigee
        LP = acosd(dot(i_vector,E_vector)/E);
        if E_vector(1,2) < 0
                LP = 360 - LP;
        end
        fprintf(fID,"Longitude of Perigee: %.4f\n", LP);
        end
    else
        if E < 0.001
            %Argument of Latitude 
            AL = acosd(dot(n_vector,R_vector)/(n*R));
            if R_vector(1,3) < 0
                AL = 360 - AL;
            end
            fprintf(fID,"Argument of Latitude: %.4f\n", AL);
        end
    end

        %Radius of Perigee 
        RP = abs(a*(1-E));
        if (RP < 6628)
            fprintf(fID,"Ballistic Missile with a Radius of Perigee: %.4f \n", RP);
        end

        if isempty(TA)
            TA = AL;
        end
        if isempty(AOP)
            AOP = 0;
        end

        shapeOfOrbit(E);
        keplersProblem(L, longitude, H, GST, E, TA, mu, a, TOF, i, RAAN, AOP)

end