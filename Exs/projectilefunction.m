% projectilefunction.m 
%
% Simulation of a projectile in a plane.
% Gravity acts down in the vertical (-y) direction.
% No energy loss (e.g., damping, drag) 
% 
% The particle is given an initial velocity (can have components in x
% and/or y directions
%
% The simulation stops when the mass touches the ground (y=0) 
%
% last modified 8/30/18  CLee
%
function projectilefunction
clear all
close all

% Define system parameters
%
g = 9.81;              % gravitational acceleration in m/s^2
mass = 10;             % mass in kg
%
% Specify magnitude and angle of initial velocity
theta    = 45;     % initial angle in degrees
v0_mag   =  8;     % magnitude of velocity in m/s

% Define state variables: 
%  z1 = x, z2 = dx/dt, z3 = y, z4 = dy/dt

% Specify initial velocity. Inital displacements are taken to be zero. 
z1_0 = 0;  
z2_0 = v0_mag*cosd(theta);% 
z3_0 = 0;    %
z4_0 = v0_mag*sind(theta);

Z_0 = [z1_0, z2_0, z3_0, z4_0];
% Define simulation parameters
t_span = [0: 0.05: 50];  % max time span for simulation 
% integration will stop at end of time span or when projectile hits ground
% (as defined by events)

options = odeset('Events', @event_stop);
[t, zout] = ode45(@projectile_fun, t_span, Z_0, options);

t(end)  %end time

%analytical solutions
y_eqn = -g/2*t.^2 + z4_0*t + z3_0;
x_eqn = z2_0*t + z1_0;

%kinetic, potential, and total menergy 
vel_mag = sqrt( zout(:,2).^2 + zout(:,4).^2 );
KE = 1/2*mass.*vel_mag.^2;
PE =  mass * g * zout(:,3);
Etotal = KE + PE;

%
% plot trajectories of upper and lower masses 
% combine phase 1 and phase 2
%
figure  %plots of time vs. hort/vert displacement
subplot(2,1,1), plot(t,zout(:,1)) 
hold
plot(t,x_eqn, 'r+')
ylabel('Horiz Displ (m)')
title('Projectile Displacements ')
subplot(2,1,2), plot (t, zout(:,3))
hold
plot(t,y_eqn, 'r+')
ylabel('Vert Displ (m)')
%
xlabel('Time (sec)')
legend ('ode45','analytical')

figure % trajectory path in space
plot(zout(:,1), zout(:,3))
hold
plot(x_eqn, y_eqn, 'r+')
title('Trajectory in Space')
ylabel('Vert Displ (m)')
xlabel('Horiz Displ (m)')
legend ('ode45','analytical')
%
% for scaling plots
% axis ([0 0.5 0 350])

figure  % KE, PE, Etotal as function of time
plot(t, KE)
hold
plot(t, PE, 'g')
plot(t, Etotal, 'r')
title('Kinetic, Potential, Total Energy')
xlabel('Time (s)')
ylabel('Energy (J)')
legend('KE', 'PE', 'total')

time_end = t(end)      % simulation time when stopped by events
horizontal_distance = zout(end,1) % horizontal distance (zout(1) = x, when stimulation stopped)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function definition in first order form
function states = projectile_fun(T, ZZ)
% solve simultaneously for x and y ( and dx/dt and dy/dt)
z1 = ZZ(2); 
z2 = 0;
z3 = ZZ(4);
z4 = -g;
%
states = [z1;z2;z3;z4];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdir] = event_stop(T,Z)
     
        % stop when Z3 = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(3);    %  ‘Events’ are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = -1;       %  Detect only events with dydt<0
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end










