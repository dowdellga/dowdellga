% AERSP 450 Project 2
% Group 6: Conor Dowdell, Gabrielle Nibert, Sebastian Valentin
% Author: Conor Dowdell

close all; clc; clear;

% Part B, Number 2- "Compute the position and velocity vectors at these
% points on the two orbits

r = 10130.5; %km (Note r_1 = r_2)
mu = 3.986E5;
a_1 = 13000;
a_2 = 7226.58;
e_1 = 0.3;
e_2 = 0.444819;
p_1 = a_1*(1-e_1^2);
f_1 = acosd( (p_1/(r*e_1) ) - (1/e_1)); % if positive, gamma is positive
p_2 = a_2*(1-e_2^2);
f_2 = acosd( (p_2/(r*e_2) ) - (1/e_2));
h_1 = sqrt(mu*p_1);
h_2 = sqrt(mu*p_2);

v_1 = [ ((mu/h_1)*e_1*sind(f_1)); (mu/h_1)*(1+e_1*sind(f_1)); 0 ]; % r, Θ, h
v_2 = [ ((mu/h_2)*e_2*sind(f_2)); (mu/h_2)*(1+e_2*sind(f_2)); 0 ];

%---------------------------------------------------------------------------

%Part C- Compute the delta-v at the intersection point to transfer from one orbit to another

gamma_1 = atan( e_1*sind(f_1) / (1+e_1*cosd(f_1)) );
gamma_2 = atan( e_2*sind(f_2) / (1+e_2*cosd(f_2)) );
delta_gamma = gamma_2 - gamma_1;

% Since there is no plane-change... inclination = 0


delta_v = norm(v_1)^2 + norm(v_2)^2 - 2*norm(v_1)*norm(v_2)*cos(delta_gamma);

fprintf('The change in velocity is equal to %.3f km/s',delta_v);

%%
% Numerical Part

close all; clc; clear;
load('IODMeasurements2.mat');

mu = 3.986E5;


a_matrix = zeros(24,1);
e_matrix = zeros(24,1);
f_matrix = zeros(24,1);
v_matrix = zeros(24,1);
r1r2r3 = zeros(24,3);

% Note:
% Asimuth & Elevation values expressed in Topocentric Equatorial Frame
% telescope is located at 31.9466 N (latitude) and 108.8977 W (Longitude)
% (31.9466 N, 108.8977 W)

% Algorithm-
% Step 1: Determine initial time, L, and R

% Test for first 3 observations!

%set = 22; % 1, 4, 7, 10, 13, 16, 19, 22 This controls the set of measurements (i.e. 4 = 2nd set, 6 equals 3rd set)

% 4 & 7 are same orbit
% 10 & 13 are same orbit

for set = 1:3:24


for i = set:1:set+2 % This controls the set of measurements (i.e. 4:1:6 is the second set, 7:1:9 is 3rd, 16:1:18 is 6th)

%set = 1;
S_i = ELEVATION(i); % degrees
A_i = AZIMUTH(i); % degrees
L_i = [cosd(A_i)*cosd(S_i);
       sind(A_i)*cosd(S_i);
       sind(S_i)];


    if i == set % 1, 4, 7, 10, 13, 16, 19, 21
    L_1 = L_i;
    R_1 = RSITES(i,:);

    % Calculate time intervals (see above code for tau_1, tau_3, etc.)
    tau_1 = TIMES(set)-TIMES(set+1); 
    tau_3 = TIMES(set+2)-TIMES(set+1);
    tau_13 = tau_3 - tau_1;

    end

    if i == set + 1
    L_2 = L_i;
    R_2 = RSITES(i,:);
    end

    if i == set + 2
    L_3 = L_i;
    R_3 = RSITES(i,:);

    end


end




%Computing 10 triple products

D_0 = dot(L_1,cross(L_2,L_3));

D_11 = dot(R_1, cross(L_2,L_3)); % 1x3 dotted with a 3x1
D_12 = dot(R_1, cross(L_1,L_3));
D_13 = dot(R_1, cross(L_1,L_2));

D_21 = dot(R_2, cross(L_2,L_3));
D_22 = dot(R_2, cross(L_1,L_3)); 
D_23 = dot(R_2, cross(L_1,L_2));

D_31 = dot(R_3, cross(L_2,L_3));
D_32 = dot(R_3, cross(L_1,L_3)); 
D_33 = dot(R_3, cross(L_1,L_2));




% Calculating Coefficients A & B

A = (1/D_0)* ( (-tau_3/tau_13)*D_12 + D_22 + (tau_1/tau_13)*D_32 );

B = (1/(6*D_0))* (-((tau_13^2)-(tau_3^2))*(tau_3/tau_13)*D_12 + ((tau_13^2)-(tau_1^2))*(tau_1/tau_13)*D_32);


a = (-A^2) - 2*A*dot(L_2,R_2)-(norm(R_2)^2);
b = -2*mu*B*(A + dot(L_2,R_2));
c = -(mu^2)*B^2;


p = [1 0 a 0 0 b 0 0 c];
r_2_array = roots(p); % 2 real (1 (+), 1 (-), 6 imaginary (2(+), 4(-))
% roots should all have same inclination and RAAN

if set == 1 || set == 13 || set == 16 || set == 19
        r_2 = r_2_array(4); % For 1st orbit & sixth orbit
elseif set == 4 ||set == 7 || set == 10
        r_2 = r_2_array(1); % this is for 2nd & 3rd & 4th orbits
elseif set == 22
        r_2 = r_2_array(3);
end

num_1 = 6.*(D_31.*(tau_1/tau_3) + D_21.*(tau_13/tau_3)).*r_2^3 + (mu.*D_31).*(tau_13^2-tau_1^2).*(tau_1/tau_3);
denom_1 = 6*r_2^3 + mu*(tau_13^2 - tau_3^2);
rho_1 = (1/D_0)*(num_1/denom_1 - D_11);

rho_2 = A + mu*B*r_2^-3;

num_3 = 6.*( D_13.*(tau_3/tau_1) - D_23.*(tau_13/tau_1)).*r_2^3 + (mu.*D_13).*(tau_13^2-tau_3^2).*(tau_3/tau_1);
denom_3 = 6*r_2^3 + mu*(tau_13^2 - tau_1^2);
rho_3 = (1/D_0)*(num_3/denom_3 - D_33);


position = [R_1' + rho_1*L_1; R_2' + rho_2*L_2; R_3' + rho_3*L_3];
r1_HGibbs = [position(1); position(2); position(3)];
r2_HGibbs = [position(4); position(5); position(6)];
r3_HGibbs = [position(7); position(8); position(9)];


coplanar_check = dot(r1_HGibbs, cross(r2_HGibbs,r3_HGibbs));
r1_mag = norm(r1_HGibbs);
r2_mag = norm(r2_HGibbs);
r3_mag = norm(r3_HGibbs);


% After computing the position vectors at t_1, t_2, & t_3, we now use
% Herrick-Gibbs method (because there are short times between observations)

del_t31 = TIMES(set+2) - TIMES(set);
del_t32 = TIMES(set+2) - TIMES(set+1);
del_t21 = TIMES(set+1) - TIMES(set);


v_2 = -del_t32*( (1/(del_t21*del_t31)) + mu/(12*r1_mag^3))*r1_HGibbs + (del_t32 - del_t21)*( (1/(del_t21*del_t32)) ...
    + mu/(12*r2_mag^3))*r2_HGibbs + del_t21*( 1/(del_t32*del_t31) + mu/(12*r3_mag^3))*r3_HGibbs;

v2_mag = norm(v_2);

v_matrix(i) = v2_mag;

fprintf('\nBelow are the orbital elements from Observations %.1f to %.1f:',set,set+2);

eps = ((dot(v_2,v_2))/2) - (mu/r2_mag); % epsilon
fprintf("\nThe Specific Energy is: %.2f km^2/s^2 ",eps);


a = -mu/(2*(eps)); % solving for the semi-major axis
fprintf("\nThe Semi Major Axis, a = %.2f km ",a);

% How can I index a???

% for example: I want to store all the different a-values
% a= [6999.93, 1; 10964.91, 2; ...]
% 1, 4, 7, 10

a_matrix(i) = a;


h = cross(r2_HGibbs,v_2); % calculating angular momentum LU^2/TU
h_mag = norm(h); % magnitude of angular momentum
e = cross(v_2,h)/mu - (r2_HGibbs/r2_mag);
e_mag = norm(e);
fprintf("\nThe eccentricity value, e = %.2f (dimensionless) ",e_mag);
e_matrix(i) = e_mag;


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
r_0 = r2_mag; % r_0 equals magnitude of r_E

f = acos( (1/e_mag) * ((p/r_0)-1) ); % true anomaly in degrees
f_deg = rad2deg(f);

f_matrix(i) = f_deg;

fprintf("\n Flight path angle, f = %.2f radians or %.2f degrees \r",f,f_deg);

fprintf('The position 2 (magnitude) = %.3f [km]\n',r2_mag);
fprintf('The velocity vector (Gibbs Method) of position 2 = [%.3f %.3f %.3f] [km/s] \n',v_2(1), v_2(2), v_2(3));

fprintf('Position vector (x,y,z) for time 1 = [%.3f %.3f %.3f] [km]\n',position(1), position(2), position(3));
fprintf('Position vector (x,y,z) for time 2 = [%.3f %.3f %.3f] [km]\n',position(4), position(5), position(6));
fprintf('Position vector (x,y,z) for time 3 = [%.3f %.3f %.3f] [km]\n',position(7), position(8), position(9));
disp('Here are the calculated roots');
disp(r_2_array);
fprintf('We chose r_2 = %.3f because it is real, positive, and because I & ω remain the same (no plane change)\n',r_2);

T_p = 2*pi*sqrt(a^3/mu);
dt = linspace(0,T_p,1000);
x0 = [r2_HGibbs, v_2];
options = odeset('reltol',1e-12,'abstol',1e-12);
[t,position_vec] = ode45(@(t,g) TwoBP(t,g,mu),dt,x0,options); % calling ODE45


figure (1)
% Creating/Plotting Spherical Earth
r_EARTH = 6378.14; % Radius of Earth [km]
[xEarth, yEarth, zEarth] = sphere(30);
surf(r_EARTH*xEarth,r_EARTH*yEarth,r_EARTH*zEarth, 'FaceColor', [0 0 1]);
%
plot3(position_vec(:,1),position_vec(:,2),position_vec(:,3),'LineWidth',3);
xlabel("x [km]");
ylabel("y [km]");
zlabel("z [km]");
title("3D Plot of 5 Orbits");
legend('Orbit 1')
grid on;
hold on;
%pause(1);

end
saveas(gcf,'Orbit_Plot.jpg', 'jpg');


%Part 7 Compute the delta-v at the intersection point to transfer from one orbit to another


k = 21; % used to compute delta v's

e_1 = e_matrix(k);
f_1 = f_matrix(k);
v_1 = v_matrix(k);
e_2 = e_matrix(k+3);
f_2 = f_matrix(k+3);
v_2 = v_matrix(k+3);


gamma_1 = atan( e_1*sind(f_1) / (1+e_1*cosd(f_1)) );
gamma_2 = atan( e_2*sind(f_2) / (1+e_2*cosd(f_2)) );
delta_gamma = gamma_2 - gamma_1;

% Since there is no plane-change... inclination = 0


delta_v = norm(v_1)^2 + norm(v_2)^2 - 2*norm(v_1)*norm(v_2)*cos(delta_gamma);

fprintf('\n Part 7- The change in velocity is equal to %.3f km/s',delta_v);





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


