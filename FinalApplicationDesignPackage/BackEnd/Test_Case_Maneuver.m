% MAE 4510_Senior_Design_Lockheed_COOLR
% Lynnane George
% 3/28/23


% This program finds the most efficient maneuver to get from one longitude 
% to another for a geostationary satellte
% Given input position and velocity vectors in the J2000 frame

clc, clear;

fileID=fopen(['D:\School\Sr Design\Test_Case_Maneuver.txt'],'w');
% Print some nice formatting for the output
fprintf(fileID,'********************************CASE_1*********************************** \n');
fprintf(fileID,'This program takes as inputs Position and Velocity Vectors and returns the COEs.  It will also find the COEs after a Time of Flight and provide siting information.  \n');
fprintf(fileID,'                       Lockheed COOLR test case\n');
fprintf(fileID,'                        by Lynnane George \n');
fprintf(fileID,'                             MAE 4510 \n');
fprintf(fileID,'                              3/5/2023 \n');


%***********************************INPUTS**********************************
% Enter satellite data in terms of
% Year 
% Month
% Day
% Time in Hrs
% Time in Mins
% Time in Secs

Year = 2022;
Month = 10.0;
Day = 31.0;
Hour = 23.0;
Min = 59.0;
Sec = 42.0;

Date=[Year Month Day Hour Min Sec];
% Input satellite Position and Velocity Vectors

x = (2.6018527189493418e+07)/1000; % km
y = -(3.1172463558978092e+07)/1000; % km
z = -(1.1367327082502998e+07)/1000; % km
xdot = (2.3565352869812309e+03)/1000 %km/s
ydot = (1.9747781263341874e+03)/1000 %km/s
zdot = -(2.5741505025935862e+01)/1000 %km/s

% Input final u
uf = 360;


% Convert Date to Julian Data
jDate = jdFromDate(Date);

% Convert J2000 to TOD
tmatrix = J2000toTOD(jDate);

% Find R and V vectors
R = tmatrix*[x;y;z];
omegaearth = [0; 0; 15*pi/(180*3600)];
V = tmatrix*[xdot;ydot;zdot];

% Print Initial Satellite Tracking Information
fprintf(fileID,'Initial Date:\n');
fprintf(fileID,'Year = %i\n',Year);
fprintf(fileID,'Month = %f\n',Month);
fprintf(fileID,'Day = %f\n',Day);
fprintf(fileID,'Hour = %f\n',Hour);
fprintf(fileID,'Minute = %f/n.',Min);
fprintf(fileID,'Seconds = %f\n',Sec);
fprintf(fileID,'Julian Date = %f\n\n', jDate);

% Define constants
% gravitational parameter for the earth
mu = 398600.5; %km^3/s^2
% epsilon = small number to check for special cases
eps = .001;

% Print Position and Velocity Vectors to file
fprintf(fileID,'R = %f  I  +  %f  J  +  %f  K\n',R(1),R(2),R(3));
fprintf(fileID,'V = %f  I  +  %f  J  +  %f  K\n\n',V(1),V(2),V(3));

% First find specific angular momentum and ascending node vectors
h=cross(R,V);
nhat=cross([0 0 1],h);

% Find energy
energy = norm(V)^2/2-mu/norm(R);

% Find the eccentricity vector
evec = ((norm(V)^2-mu/norm(R))*R-dot(R,V)*V)/mu;
e = norm(evec)

% a = semi-major axis
if abs(e-1.0)>eps
   a = -mu/(2*energy);
   p = a*(1-e^2);
else
   % this is a parabola
   p = norm(h)^2/mu;
   a = inf;
end

% print first two COEs to output file
fprintf(fileID,'The COEs for the spacecraft are:  \n');
fprintf(fileID,'a:  semi-major axis =  %f  \n',a);
fprintf(fileID,'e:  eccentricity =  %f  \n',e);

% find inclination
i = acos(h(3)/norm(h))*180/pi;
fprintf(fileID,'i:  inclination =  %f   \n',i);

% check for non-equatorial case
if i>eps

% Calculate RAAN
Omega = acos(nhat(1)/norm(nhat))*180/pi;

% Perform half plane check
if nhat(2)<0
   Omega = 360-Omega;
end

fprintf(fileID,'Omega:  RAAN =  %f     \n', Omega);

% Check for elliptical case
if e > eps
    
    % calculate argument of perigee
    argp = acos(dot(nhat,evec)/(norm(nhat)*e))*180/pi;
      % perform half-plane check
      if evec(3)<0
    argp = 360-argp;
      end
      
    % calculate true anomaly
    nu = acos(dot(evec,R)/(e*norm(R)))*180/pi;
    % perform half-plane check
    if dot(R,V)<0
       nu = 360 - nu;
    end

    fprintf(fileID,'Argp:  argument of perigee =  %f     \n',argp);
    fprintf(fileID, 'Nu:  true anomaly  =  %f        \n',nu);
    
  else
    % this is a circular orbit
    % calculate alternate COE u = argument of latitude
    % measure from n to R
    u = acos(dot(nhat,R)/(norm(nhat)*norm(R)))*180/pi;
      % perform half-plane check
      if (R(3)<0)
        u = 360 - u;
      end
    fprintf(fileID, 'u:  argument of latitude = %f    \n',u);  
  end
  
else
    if (e>eps)    
      % calculate alternate COE longper =  longitude of perigee
         longper = acos(evec(1)/e)*180/pi;
      % perform half-plane check
      if(evec(2)<0)
        longper = 360 - longper;
      end
    fprintf(fileID,'longper = longitude of perigee = %f     \n',longper);
    
    else
        % this is a circular equatorial case
        % calculate true longitude
        % measure from I to R
        truelong=acos(R(1)/norm(R))*180/pi;
        % perform half-plane check
        if (R(2)<0) 
            truelong=360-truelong;
        end
        fprintf(fileID,'l:  true longitude = %f      \n',truelong);
     end
     
end

% print appropriate trajectory type
if abs(e-1.0)<eps
    fprintf(fileID,'This is a parabolic orbit\n');
   fprintf(fileID,'\n*******************************************************************');
else 
    if abs(e)<eps
         fprintf(fileID,'This is a circular orbit\n');
        fprintf(fileID,'\n*******************************************************************');
    
    else if abs(e)>1
            fprintf(fileID,'This is a hyperbolic orbit\n');
            fprintf(fileID,'*******************************************************************');
        else
            fprintf(fileID,'This is an elliptical orbit\n');
            fprintf(fileID,'*******************************************************************');
        end
end
end

if abs(i)<eps
        fprintf(fileID,'This is an equatorial orbit \n\');
        fprintf(fileID,'*******************************************************************');
end


% Maneuver
for n = 1:10
    e(n) = n/10;
    vi (n) = acos((a*(1-e(n)^2)/(a)-1));
    Ra(n) = a*(1+e(n));
    at = (a + Ra(n))/2;
    Et = -mu/(2*at);
    Vt1 = sqrt(2*(mu/norm(R)-Et));
    phi = atan(e(n)*sin(vi(n)/(1+e(n)*cos(vi(n)))));
    Dv(n) = sqrt(Vt1^2+(norm(V))^2-2*Vt1*norm(V)*cos(phi));

 end
 

fclose(fileID);

