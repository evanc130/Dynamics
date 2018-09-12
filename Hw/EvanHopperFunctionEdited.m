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
function hopperfunction
clear all
close all

% Define system parameters
%
g = 9.81;              % gravitational acceleration in m/s^2
ml = .010;             % mass in kg 10GRAMS
mu = .020;             % 20 GRAMS
k = 9810;              %N/m^2
L0=.050;               %m
d=.005;                %m
%

% Define state variables: 
%  z1 = x, z2 = dx/dt, z3 = y, z4 = dy/dt

% Specify initial velocity. Inital displacements are taken to be zero. 
z1_0 = L0-d;  
z2_0 = 0;		% vu
z3_0 = 0;    					%
z4_0 = 0;		% vl

Z_01 = [z1_0, z2_0, z3_0, z4_0];

% Define simulation parameters
T_span1 = [0: 0.001: 50];
 

options1 = odeset('Events', @event_switch);
options2 = odeset('Events', @event_stop);
[t1, zout1] = ode45(@hopper_fun1, T_span1, Z_01, options1);

%Start the figure
figure 
plot(t1,zout1(:,1));
hold on;
plot(t1,zout1(:,3));

Z_02 = [zout1(end,1), zout1(end,2), zout1(end,3), zout1(end,4)];

T_span2 = [t1(end): 0.001: 50];
   
[t2,zout2] = ode45(@hopper_fun2,T_span2, Z_02, options2);

t1(end)+t2(end)  %end time
t2(end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %analytical solutions
   

% 
% %kinetic, potential, and total energy 
% vel_mag = sqrt( zout1(:,2).^2 + zout1(:,4).^2 );
% KE = 1/2*mass.*vel_mag.^2;
% PE =  mass * g * zout1(:,3);
% Etotal = KE + PE;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% print time values
% plot(t1,zout1(1));
% hold on;
% plot(t2,zout2(1));
% plot(t1,zout1(3));
% plot(t2,zout2(3));

% trajectory path in x-y space
plot(t1, zout1(:,1));
hold on;
plot(t1, zout1(:,3));
hold on;
plot(t2,zout2(:,1));
hold on;
plot(t2,zout2(:,3));
% title('Trajectory in Space')
% ylabel('Vert Displ (m)')
% xlabel('Horiz Displ (m)')
% legend ('ode45')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function: equations of motion in first order form
function dzdt1 = hopper_fun1(T,Z)
% solve simultaneously for x and y,  and dx/dt and dy/dt
% Z(1)= z1=x, Z(2)= z2=dx/dt, Z(3)= z3=y, Z(4)= z4=dy/dt
dz1dt = Z(2); 
dz2dt = -g+k/mu*(L0-Z(1));
dz3dt = Z(4);
dz4dt = 0;
%
dzdt1 = [dz1dt;dz2dt;dz3dt;dz4dt];
%
end

function dzdt2 = hopper_fun2(T,Z)
% solve simultaneously for x and y,  and dx/dt and dy/dt
% Z(1)= z1=x, Z(2)= z2=dx/dt, Z(3)= z3=y, Z(4)= z4=dy/dt
dz1dt2 = Z(2); 
dz2dt2 = -g-k/mu*((Z(1)-Z(3))-L0);
dz3dt2 = Z(4);
dz4dt2 = -g+k/ml*((Z(1)-Z(3))-L0);
%
dzdt2 = [dz1dt2;dz2dt2;dz3dt2;dz4dt2];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdir] = event_switch(T,Z)
     
        % stop when Z(3)= z3= y = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(1)-Z(3)-L0;    %  ‘Events’ are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = 1;       %  Detect only events with dydt<0
end
 
function [eventvalue,stopthecalc,eventdir] = event_stop(T,Z)
     
        % stop when Z(3)= z3= y = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(3);    %  ‘Events’ are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = 0;       %  Detect only events with dydt<0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end










