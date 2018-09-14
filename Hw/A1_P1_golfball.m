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
function A1_P1_golfball
clear all
close all

% Define system parameters
%
g = 9.81;              % gravitational acceleration in m/s^2
m = 0.0459;             % mass in kg
rho = 1.29;               %density in kg/m^3
Cd = 0.25;              %Drag Coefficient
% Cd = 0;                 %Model Validation: disable wind
D = 4.27/100;           %Diameter in m

c=pi*rho*Cd*D^2/8/m;
%
% Specify magnitude and angle of initial velocity
x_0=0;
y_0=0;
z_0=0;

metric=2.2;

v=130/metric
th=45;

vx_0=v*cosd(th);
vy_0=0;
vz_0=v*sind(th);

vx_w=0/metric;
vy_w=15/metric;
vz_w=0/metric;

% for i=1:3
% if i==1
%         vx_w=15/metric;
%         vy_w=0/metric;
%         vz_w=0/metric;
% else if i==2
%         vx_w=-15/metric;
%         vy_w=0/metric;
%         vz_w=0/metric;
%     else if i==3
%         vx_w=0/metric;
%         vy_w=15/metric;
%         vz_w=0/metric;
%         end
%     end
% end
%         
% 
% 
% 
% Z_0 = [x_0, vx_0, y_0, vy_0, z_0, vz_0];
% % Define simulation parameters
% t_span = [0: 0.05: 50];  % max time span for simulation 
% % integration will stop at end of time span or when projectile hits ground
% % (as defined by events)
% 
% options = odeset('Events', @event_stop);
% [t, zout] = ode45(@projectile_fun, t_span, Z_0, options);
% 
% t(end)  %end time
% 
% figure(1);
% plot3(zout(:,1),zout(:,3),zout(:,5));
% axis equal;
% title('Ball Path');
% xlabel('x(i)');
% ylabel('y(j)');
% zlabel('z(k)');
% hold on;
% % end
% legend('Tailwind','Headwind','Crosswind');



Z_0 = [x_0, vx_0, y_0, vy_0, z_0, vz_0];
% Define simulation parameters
t_span = [0: 0.05: 50];  % max time span for simulation 
% integration will stop at end of time span or when projectile hits ground
% (as defined by events)

options = odeset('Events', @event_stop);
[t, zout] = ode45(@projectile_fun, t_span, Z_0, options);

t(end)  %end time

figure(1);
plot3(zout(:,1),zout(:,3),zout(:,5));
axis equal;
title('Ball Path');
xlabel('x(i)');
ylabel('y(j)');
zlabel('z(k)');

a=zout(end,1)

figure(2);
subplot(3,1,1), plot (t, zout(:,1)), hold on;
subplot(3,1,1), plot(t, zout(:,3)), hold on;
subplot(3,1,1), plot(t,zout(:,5)), hold on;
legend('x','y','z')
title('Position');
xlabel('time(s)');
ylabel('Position (m)');

subplot(3,1,2), plot (t, zout(:,2)), hold on;
subplot(3,1,2), plot(t, zout(:,4)), hold on;
subplot(3,1,2), plot(t,zout(:,6)), hold on;
legend('x','y','z')
title('Velocity');
xlabel('time(s)');
ylabel('velocity (m/s)');

%If I can cheat and use ode45, then I can cheat and derrive acceleration numerically
ax=diff(zout(:,2));
ay=diff(zout(:,4));
az=diff(zout(:,6));

subplot(3,1,3), plot (t(2:end), ax), hold on;
subplot(3,1,3), plot(t(2:end), ay), hold on;
subplot(3,1,3), plot(t(2:end), az), hold on;
legend('x','y','z')
title('Acceleration');
xlabel('time(s)');
ylabel('Acceleration (m/s^2)');










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function definition in first order form
function states = projectile_fun(T, ZZ)
% solve simultaneously for x and y ( and dx/dt and dy/dt)

V=sqrt((ZZ(2)-vx_w)^2+(ZZ(4)-vy_w)^2+(ZZ(6)-vx_w)^2);

z1 = ZZ(2); 
z2 = -c*V*(ZZ(2)-vx_w);
z3 = ZZ(4);
z4 = -c*V*(ZZ(4)-vy_w);
z5 = ZZ(6);
z6 = -g-c*V*(ZZ(6)-vz_w);
%
states = [z1;z2;z3;z4;z5;z6];
%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eventvalue,stopthecalc,eventdir] = event_stop(T,Z)
     
        % stop when Z3 = 0 (mass hits the ground in y-dir)
        eventvalue  =  Z(5);    %  ‘Events’ are detected when eventvalue=0
        stopthecalc =  1;       %  Stop if event occurs
        eventdir    = -1;       %  Detect only events with dydt<0
end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end










