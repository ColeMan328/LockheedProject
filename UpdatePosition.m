% User inputs
tspan = input('Enter time span in seconds: ');
r0 = [26018.5271, -31172.4636, -11367.327082];
v0 = [2.3565, 1.97477, -0.0257];

% Constants
mu = 398600.5; % Gravitational parameter for Earth

% Define the differential equation
f = @(t,Y) [Y(4); Y(5); Y(6); ...
            -mu/norm(Y(1:3))^3*Y(1); ...
            -mu/norm(Y(1:3))^3*Y(2); ...
            -mu/norm(Y(1:3))^3*Y(3)];

% Combine the initial position and velocity vectors
Y0 = [r0; v0];

% Solve the differential equation using ode45
[t,Y] = ode45(f, [0 tspan], Y0);

% Extract the position and velocity vectors from the solution
r = Y(:,1:3);
v = Y(:,4:6);

% Plot the trajectory
figure;
plot3(r(:,1),r(:,2),r(:,3),'linewidth',2);
grid on;
title('Satellite Trajectory');
xlabel('x');
ylabel('y');
zlabel('z');

% Display the final position and velocity vectors
fprintf('Final position vector: [%f %f %f]\n', r(end,:)); 
fprintf('Final velocity vector: [%f %f %f]\n', v(end,:));