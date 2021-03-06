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
ml = 0.01;             % mass in kg
mu = 0.02;
k = 9810;               %N/m
L0=.05;                %mm
d=.005;                    %mm
Tstep=0.0001;
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
T_span1 = [0: Tstep: 5];
% time range for simulation with specified time step
% integration will stop at end of time span or when projectile hits ground
% (as defined by events)

options1 = odeset('Events', @event_switch);
options2 = odeset('Events', @event_stop);

[t1, zout1] = ode45(@hopper_fun1, T_span1, Z_01, options1);

Z_02 = [zout1(end,1), zout1(end,2), zout1(end,3), zout1(end,4)];

T_span2 = [t1(end): Tstep: 1];

[t2,zout2] = ode45(@hopper_fun2,T_span2, Z_02, options2);

  %end time
t2(end)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal force
n1=-(-g*mu-g*ml-k*(L0-zout1(:,1)));
% n2=-(-g/mu-g/ml+k/mu/ml*(L0-zout2(:,1)-zout2(:,3)));

figure (1)
plot(t1,n1);
% hold on;
% plot(t2,n2);
title('Normal Force');  xlabel('time (s)'); ylabel('Force (N)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kinetic, potential, and total energy
v1u=zout1(:,2); v1l=zout1(:,4); v2u=zout2(:,2); v2l=zout2(:,4);

%Phase 1
KE1u=1/2*(mu).*v1u.^2;
KE1l=1/2*(ml).*v1l.^2;
KE1 = KE1u+KE1l;

PE1u =  mu * g * zout1(:,1);
PE1l =  ml * g * zout1(:,3);
PE1 = PE1u + PE1l;



%Phase 2
KE2u=1/2*(mu).*v2u.^2;
KE2l=1/2*(ml).*v2l.^2;
KE2 = KE2u+KE2l;

PE2u =  mu * g * zout2(:,1);
PE2l =  ml * g * zout2(:,3);
PE2 = PE2u + PE2l;

%spring
Espring1=abs(1/2*k*(zout1(:,1)-zout1(:,3)-L0).^2);
Espring2=abs(1/2*k*(zout2(:,1)-zout2(:,3)-L0).^2);

Etotal1 = KE1 + PE1+Espring1;
Etotal2 = KE2 + PE2+Espring2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)

plot(t1,Etotal1);
hold on;
plot(t1,KE1);
hold on;
plot(t1,PE1);
hold on;
plot(t1,Espring1);
title('Energy Phase 1');
legend('total','kinetic','potential','spring');


figure (3)

plot(t2,Etotal2);
hold on;
plot(t2,KE2);
hold on;
plot(t2,PE2);
hold on;
plot(t2,Espring2);
title('Energy Phase 2');
legend('total','kinetic','potential','spring');



%%%%%%%%
% trajectory path in x-y space
figure (4) 
plot(t1, zout1(:,1));
hold on;
plot(t1, zout1(:,3));
hold on;
plot(t2,zout2(:,1));
hold on;
plot(t2,zout2(:,3));
xlabel('time (s)'); ylabel('Height (m)');
title('Hopper Top Heavy');
title('Hopper Bottom Heavy');

title('Trajectory in Space')
ylabel('Vert Displ (m)')
xlabel('Horiz Displ (m)')
legend ('ode45')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function: equations of motion in first order form
function dzdt1 = hopper_fun1(T,Z)
% solve simultaneously for x and y,  and dx/dt and dy/dt
% Z(1)= z1=yu, Z(2)= z2=dyu/dt, Z(3)= z3=yl, Z(4)= z4=dyl/dt
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
% Z(1)= z1=yu, Z(2)= z2=dyu/dt, Z(3)= z3=yl, Z(4)= z4=dyl/dt
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
        eventvalue  =  -(-g/mu-g/ml-k/mu/ml*(L0-Z(1)));    %  �Events� are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = 0;       %  Detect only events with dydt<0
end
 
function [eventvalue,stopthecalc,eventdir] = event_stop(T,Z)
     
        % stop when Z(3)= z3= y = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(3);    %  �Events� are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = -1;       %  Detect only events with dydt<0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end










