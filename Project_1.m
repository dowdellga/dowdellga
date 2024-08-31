%
% AERSP 450 Project 1
% Conor Dowdell, Gabrielle Nibert, Sebastian Valentin

clear; clc; close all;


%syms GMST

% defining knowns
a = 6652;
e_mag = 0.2;
h = 50452.2;
mu = 3.986*10^5;
t = 1.5; % 5400 secs
r_p = 5321.6; % radius at periapsis
phi = 45; % degreees North (min latitude is 41 degrees)
lambda = -78; % west
cap_omg = 60;
GMST = 230.1; % found via analytical solution
I = (45); % chosen because c(I) & s(I) both = sqrt(2)/ 2
del = 1e-8;
r0 = [4759.431; 31.56621; 4759.536]';
v0 = [1.3825; 7.9809; 3.0614]';
r0_mag = norm(r0);
v0_mag = norm(v0);
f = 0;
E_0 = 2*atan(sqrt((1-e_mag)/(1+e_mag))*tan(f/2));

% Building Earth Centered Earth-Fixed frame
x = r_p*cosd(phi)*cosd(lambda);
y_ECEF = r_p*cosd(phi)*sind(lambda);
z_ECEF = r_p*sind(phi);
r_ECEF = vpa([x; y_ECEF; z_ECEF],7); % r vector in the ECEF frame
ECEF_2_ECI = vpa([cosd(GMST) -sind(GMST) 0; sind(GMST) cosd(GMST) 0; 0 0 1]); % omg is in radians in this form
r_ECI = vpa((ECEF_2_ECI * r_ECEF),7); % vpa,# where '#' represents num of sig figs


% 4 (b)- Also write the position vector of the spacecraft in the Orbital plane

r_orb = [r_p, 0, 0]'; % r/hat, theta/hat, h/hat

% defining Rotations from Ground tracks handout

R_3_cap_omg = [cosd(cap_omg) sind(cap_omg) 0;
       -sind(cap_omg) cosd(cap_omg) 0;
           0            0         1];
R_1 = [ 1       0            0;
      0 cosd(I) sind(I);
      0 -sind(I) cosd(I)];
R_3_omg = [cosd(GMST) sind(GMST) 0;
       -sind(GMST) cosd(GMST) 0;
           0            0         1];
C_NO = R_3_cap_omg*R_1*R_3_omg;

C_ON = vpa(C_NO',4);

% after solving for the orbital velocity using the initial conditions I
% know...

v_orb = [0 9.481 0]; % this is the v_orb I calculated analytically

v_ECI = vpa(C_ON * v_orb',5); % from the "Transformation from inertial position/velocity handout"
v = norm(v_ECI);

v_ECEF = ECEF_2_ECI' * v_ECI;

v_ECEF = v_ECEF';

% RV2OE(r,v,mu); % calls RV2OE func

%ECI2ECEF(r_ECI,v_ECI,cap_omg, t, GMST, mu); %calls ECI2ECEF func


% 5(a)
% Propagate and Plot the orbit of the satellite in the ECI frame. Propagate using ODE45.
T = 24*60*60; % 1.5 hours *3600 (sec/hr)
mu = 3.986*10^5;
r_ECI = [835; -3929; 3491]';
v_ECI = [1.38; 7.98; 3.06]';
r_ECEF = [989.5631 -4655.528 4759.536];
v_ECEF = [-8.2107, 3.3520, 3.3520];

dt = linspace(0,5400,1000);
x0 = [r_ECI, v_ECI];
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,position_vec] = ode45(@(t,g) TwoBP(t,g,mu),dt,x0,options); % calling ODE45





figure (1)
plot3(position_vec(:,1),position_vec(:,2),position_vec(:,3),'LineWidth',3);
xlabel("x [km]");
ylabel("y [km]");
zlabel("z [km]");
title("3D Plot of Orbit in ECEF Frame");
grid on
hold on;
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(30);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceAlpha', 0.75);

figure (2)
plot(position_vec(:,1),position_vec(:,2))% x-y plane view
xlabel("x [km]");
ylabel("y [km]");
title("X-Y Plane View");
hold on;
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(30);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceColor', [0 0 1]);

figure (3)
plot(position_vec(:,2),position_vec(:,3)) % y-z view
xlabel("y [km]");
ylabel("z [km]");
title("Y-Z Plane View");
hold on;
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(30);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceColor', [0 0 1]);

figure (4)
plot(position_vec(:,1),position_vec(:,3)) % x-z view
xlabel("x [km]");
ylabel("z [km]");
title("X-Z Plane View");
hold on;
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(30);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceColor', [0 0 1]);


figure (5)
load('topo.mat','topo');
topoplot = [topo(:,181:360),topo(:,1:180)];
contour(-180:179,-90:89,topoplot,[0,0],'black','linewidth',1);
grid on
grid minor
axis equal
hold on
x = position_vec(:,1);
y = position_vec(:,2);
z = position_vec(:,3);
%long = -78;
%lat = 41;
long = atan2d(y, x);
lat = atan2d( z, ( sqrt(x.^2 + y.^2) ) );

plot(long,lat,'m','LineWidth',2)
% <- ENTER THE LONGITUDE AND LATITUDE YOU HAVE COMPUTED HERE


xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title('Satellite Ground Track');
set(gca,'FontSize',18)


wdmap = load('coast');
wdlon = wdmap.coastlon;
wdlat = wdmap.coastlat;
set(gca,'XLim',[-180 180],'YLim',[-90 90], ...
    'XTick',[-180 -120 -60 0 60 120 180], ...
    'Ytick',[-90 -60 -30 0 30 60 90]);
plot(wdlon,wdlat,lambda,phi,'r.');axis([-200 200 -100 100]);
xlabel('Longitude ( \circ )'),ylabel('Latitude ( \circ )')
legend ('World Map','Satellite Track');
grid on;


% Creating Figure
figure (6)
%hold on
title('Two-Body Trajectory', 'Interpreter', 'Latex')
xlabel('x', 'Interpreter', 'Latex')
ylabel('y', 'Interpreter', 'Latex')
zlabel('z', 'Interpreter', 'Latex')
axis equal
grid minor
view(30, 30)
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(25);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceColor', [0 0 1]);
% Plotting Trajectory
plot3(position_vec(:,1), position_vec(:,2), position_vec(:,3), 'k')
hold off




delta_t = t(2)-t(1);
E_new = kepler(a,e_mag,mu,delta_t,del);


for i = 50
   
    % calc M for current time interval
       M = sqrt(mu/a^3)*(delta_t);
    % call newton raphson for each time int
      E_new = kepler(a,e_mag,mu,delta_t,del);

    % calc F & G with updated magnitudes
    F = (1- (a/r0_mag)*(1-cos(E_new-E_0)))*r0;
    G = (delta_t + ( ((sin(E_new-E_0)-(E_new-E_0)) / h)))*v0;
    
    % Calc r_1 from F&G
    r_1 = F.*r0 + G.*v0;


    % find mag from r_1

    % calc F_dot & G_dot from updated mag.

    %calc v_1 new

    % set r_0 and v_0 to r_1 & v_1

end


% Main Code - Problem 1
%--- Enter code here --

% Main Code - Problem 2
%--- Enter code here --

%------------------------------
% FUNCTIONS!!
%------------------------------




% Function 1: RV2OE

function [a,e,I,RAAN,AOP,f] = RV2OE(r,v,mu)

%fprintf("\n Given %.3f [km] & %.3f [km/s]...",r,v)
phi = 45;
lambda = 78;
GMST = 282.38;
cap_omg = 60;
I = 45;
x_ECF = r*cosd(phi)*cosd(lambda);
y_ECF = r*cosd(phi)*sind(lambda);
z_ECF = r*sind(phi);
ECEF = vpa([x_ECF; y_ECF; z_ECF],7);
ECEF_2_ECI = vpa([cosd(GMST) -sind(GMST) 0; sind(GMST) cosd(GMST) 0; 0 0 1]); % omg is in radians in this form

r_ECI = ((ECEF_2_ECI * ECEF)); % vpa,# where '#' represents num of sig figs
v_orb = [0 8.659, 0];



R_3_cap_omg = [cosd(cap_omg) sind(cap_omg) 0;
       -sind(cap_omg) cosd(cap_omg) 0;
           0            0         1];
R_1 = [ 1       0            0;
      0 cosd(I) sind(I);
      0 -sind(I) cosd(I)];
R_3_omg = [cosd(GMST) sind(GMST) 0;
       -sind(GMST) cosd(GMST) 0;
           0            0         1];
C_NO = R_3_cap_omg*R_1*R_3_omg;

C_ON = vpa(C_NO',4);
v_ECI = vpa(C_ON * v_orb',5);disp('Given Position vector =');
disp(r_ECI);
disp('Velocity vector =');
disp(v_ECI);
t = 5400; % 1.5 hours

r_mag = norm(r_ECI);
v_mag = norm(v_ECI);

eps = ((dot(v_ECI,v_ECI))/2) - (mu/r_mag); % epsilon
fprintf("The Specific Energy is: %.2f km^2/s^2 ",eps);

a = -mu/(2*(eps)); % solving for the semi-major axis
fprintf("\nThe Semi Major Axis, a = %.2f km ",a);


h = cross(r_ECI,v_ECI); % calculating angular momentum LU^2/TU
h_mag = norm(h); % magnitude of angular momentum
e = cross(v_ECI,h)/mu - (r_ECI/r_mag);
e_mag = norm(e);
fprintf("\nThe eccentricity value, e = %.2f (dimensionless) ",e_mag);

I = acos((h(3))/h_mag); % don't need a quad. check
I_deg = rad2deg(I); %angle btw equatorial plane & orbital plane
fprintf("\nThe Inclination Angle is: %.2f (radians) or %.2f degrees",I,I_deg);

% Before calculating Ω & ω, we need to calc. node vector 'n'

h_k = [0, 0, h(3)]; % k component of angular momentum
kxh = cross(h_k,h); % value of h_k crossed with h
n = (kxh)/(norm(kxh)); %node vector
n_J = n(2); % J component
n_I = n(1); % I component

if n(1) < 0
    Omega = atan(n(2)/n(1)) + pi; %RAAN
else
    Omega = atan(n(2)/n(1));
end
Omega_deg = rad2deg(Omega);
fprintf("\n Ω is %.2f radians or %.2f degrees ",Omega,Omega_deg);

% lil omega

lil_omega = acos(dot(e,n)/e_mag);

if e(3) < 0
    lil_omega = -lil_omega; %quad check if k component is less than zero

end

lil_omega_deg = rad2deg(lil_omega);
fprintf("\n ω is %.2f radians or %.2f degrees ",lil_omega,lil_omega_deg);

p = a*(1-e_mag^2); % finding semi-latus rectum value
r_0 = r_mag; % r_0 equals magnitude of r_E

f = acos( (1/e_mag) * ((p/r_0)-1) ); % true anomaly in degrees
f_deg = rad2deg(f);

fprintf("\n Flight path angle, f = %.2f radians or %.2f degrees \r",f,f_deg);

end


% Function #3

function [r_ECEF,v_ECEF] = ECI2ECEF(r_ECI,v_ECI,cap_omg, t, GMST, mu)

mu = 3.986*10^5;
disp('Given Position vector in the inertial frame =');
disp(r_ECI);
disp('Velocity vector in the inertial frame =');
disp(v_ECI);
disp('Ω =');
disp(cap_omg);
fprintf("\n Time in seconds, t = %.2f hours\r",t);
fprintf("\n GMST in degrees, γ = %.2f degrees\r",GMST);

ECEF_2_ECI = vpa([cosd(GMST) -sind(GMST) 0; sind(GMST) cosd(GMST) 0; 0 0 1]); % omg is in radians in this form
ECI_2_ECEF = ECEF_2_ECI'; % transpose of a DCM = inverse
r_ECEF = ECI_2_ECEF * r_ECI;
v_ECEF = ECI_2_ECEF * v_ECI;

disp('The resultant position vector in the ECEF frame =');
disp(r_ECEF);
disp('The resultant velocity vector in the ECEF frame =');
disp(v_ECEF);

end

%{
%% Function 2: OE2RV
function [output] = function_name(input)
{
% Enter code Here
}

In addition to the function for analytic part, include the following functions:

1) FGFunc - A function that codes up the F & G Functions of Lagrange. The inputs to this function are
initial position and velocity vector, and the time elapsed (∆t), the output is the final position and velocity
vector

%}

function [r,v] = FGFunc(r0,v0,delta_t,mu)

mu = 3.986*10^5;
t = 5400;
a = 6652;
r0 = [4759.431; 31.56621; 4759.536]';
v0 = [1.3825; 7.9809; 3.0614]';
r0_mag = norm(r0);
v0_mag = norm(v0);

% same code from RV2OE
RV2OE(r0,v0);

% For Elliptic Orbits -> F&G #2

F = (1- (a/r0_mag)*(1-cos(E-E_0)))*r0;
G = (delta_t + ( ((sin(E-E_0)-(E-E0)) / h)))*v0;
F_dot = -sqrt(mu*a/(r_mag*r0))*sin(E-E0)*r0;
G_dot = (1-(a/r_mag))*(1-cos(E-E0))*v0;

E_0 = 2*atan(sqrt((1-e_mag)/(1+e_mag))*tan(f/2));
% note that theta needs to be in degrees in this equation

t_p = t - sqrt(a^3/mu)*(E_0-e_mag*sin(E_0));
fprintf("\nThe time at periapsis is: %.2f seconds",t_p);

end


function E_new = kepler(a,e,mu,delta_t,del)
    a = 6652;
    M = sqrt(mu/a^3)*(delta_t); % M = 2.388;
    E_old = sqrt((mu)/(a^3))*(delta_t);
    err = 1;
    while err > del
     E_new = E_old-((M+e*sin(E_old)-E_old)/(e*cos(E_old)-1));
     err = abs(E_new-E_old);
     E_old = E_new;
    end

end

% TwoBP - A function that codes the equations of motion of the Two Body Problem for use with ode45

function [dx] = TwoBP(t,g,mu)

% call the function by saying TwoBP(time, position vector, mu value) in the
% command window OR in the script like I did with RV2OE above

x = g(1);
y = g(2);
z = g(3);
x_dot = g(4);
y_dot = g(5);
z_dot = g(6);

x_dub_dot = (-mu*x)/( (x^2) + (y^2) + (z^2) )^1.5; % found on page 9 of orbital description handout
y_dub_dot = (-mu*y)/( (x^2) + (y^2) + (z^2) )^1.5;
z_dub_dot = (-mu*z)/( (x^2) + (y^2) + (z^2) )^1.5;

dx = [x_dot; y_dot; z_dot; x_dub_dot; y_dub_dot; z_dub_dot];

end


