% projectilefunction.m 
%
% Simulation of a projectile in a plane.
% Gravity acts down in the vertical (-y) direction.
% No energy loss (e.g., damping, drag) 
% 
% The particle is given an initial velocity (can have components in x
% and/or y directions) in terms of magnitude and direction (angle wrt horizontal)
%
% The simulation stops when the mass touches the ground (y=0) 
%
% last modified 9/2/18  CLee
%
function projectilefunction
clear all
close all

% Define system parameters
%
g = 9.81;              % gravitational acceleration in m/s^2
mass = 10;             % mass in kg
%
% Specify magnitude and angle of initial velocity vector 
theta    = 45;     % initial angle (degrees) of velocity vector wrt to horizontal
v0_mag   =  8;     % magnitude of velocity in m/s

% Define state variables: 
%  z1 = x, z2 = dx/dt, z3 = y, z4 = dy/dt

% Specify initial velocity. Inital displacements are taken to be zero. 
z1_0 = 0;  
z2_0 = v0_mag*cosd(theta);		% component of velocity in x-dir
z3_0 = 0;    					%
z4_0 = v0_mag*sind(theta);		% component of velocity in y-dir

Z_0 = [z1_0, z2_0, z3_0, z4_0];
% Define simulation parameters
T_span = [0: 0.05: 50];  % time range for simulation with specified time step
% integration will stop at end of time span or when projectile hits ground
% (as defined by events)

options = odeset('Events', @event_stop);
[t, zout] = ode45(@projectile_fun, T_span, Z_0, options);

t(end)  %end time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%analytical solutions
y_eqn = -g/2*t.^2 + z4_0*t + z3_0;
x_eqn = z2_0*t + z1_0;

%kinetic, potential, and total energy 
vel_mag = sqrt( zout(:,2).^2 + zout(:,4).^2 );
KE = 1/2*mass.*vel_mag.^2;
PE =  mass * g * zout(:,3);
Etotal = KE + PE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots
figure (1) 		%plots of time vs. hort/vert displacement
subplot(2,1,1), plot(t,zout(:,1)) 
hold
plot(t,x_eqn, 'r+')  %analytical
ylabel('Horiz Displ (m)')
title('Projectile Displacements ')
subplot(2,1,2), plot (t, zout(:,3))
hold
plot(t,y_eqn, 'r+') %analytical
ylabel('Vert Displ (m)')
%
xlabel('Time (sec)')
legend ('ode45','analytical')

figure (2) % trajectory path in x-y space
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

figure (3) % Energies: KE, PE, E-total as function of time
plot(t, KE)
hold
plot(t, PE, 'g')
plot(t, Etotal, 'r')
title('Kinetic, Potential, Total Energy')
xlabel('Time (s)')
ylabel('Energy (J)')
legend('KE', 'PE', 'total')

% print time values
time_end = t(end)      				% simulation time when stopped by events
horizontal_distance = zout(end,1) 	% horizontal distance (zout(1) = x, when stimulation stopped)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function: equations of motion in first order form
function dzdt = projectile_fun(T,Z)
% solve simultaneously for x and y,  and dx/dt and dy/dt
% Z(1)= z1=x, Z(2)= z2=dx/dt, Z(3)= z3=y, Z(4)= z4=dy/dt
dz1dt = Z(2); 
dz2dt = 0;
dz3dt = Z(4);
dz4dt = -g;
%
dzdt = [dz1dt;dz2dt;dz3dt;dz4dt];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdir] = event_stop(T,Z)
     
        % stop when Z(3)= z3= y = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(3);    %  ‘Events’ are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = -1;       %  Detect only events with dydt<0
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end










