%this scripts inputs the range and range rate vectors (IN CANONICAL UNITS)
%of a satellite in a Tundra orbit around Earth. It also inputs the location
%of the radar site, time zone of the observation, time of the observation, 
%and time of flight between two points.
%With these inputs, the script:

%1. Converts the intial range and range rate vectors to ECI position and
%velocity vectors.

%2. Calculates the orbital elements of the satellite.

%3. Confirms the satellite is in a tundra orbit.

%4. Determines the true anomaly at the new point.

%5. Calculates the ECI position and velocity at this new point.

%Anthony D'Amico
%Febuary 25th, 2022


rho = [4.37923, 2.37638, 0.025397];
rh = rho + [0,0,1];

rhoDot = [0.065168, -0.16369, -0.380347];

%gravitational parameter of Earth
mu = 3.986004 * 10^5;

%date of Febuary 10th, 2022
m = 2;
d = 10;
y = 2022;

%convert 1735 MST to UT
hour = 17;
minute = 35;
tz = 'MST';

htof = 7;
mtof = 25;

if tz == 'EST'
    hour = hour + 5;
elseif tz == 'CST'
    hour = hour + 6;
elseif tz == 'MST'
    hour = hour + 7;
elseif tz == 'PST'
    hour = hour + 8;
end
    
time = hour + (minute / 60);


%location of the radar site
longitude = 104.695; %west
latitude = 38.8236; %north
direction = 'W';

%local sidereal time
J0 = ((367 * y) - (floor(7 * (y + (floor((m + 9)/12))))/4) + floor((275 * m) / 9) + d + 1721013.5);
T0 = (J0 - 2451545) / 36525;

thetag0 = 100.4606184 + (36000.77004 * (T0)) + (0.000387933 * (T0 ^ 2)) + - ((2.583 * 10 ^ -8) * (T0 ^ 3));
step1 = floor(thetag0 / 360);
thetag0 = thetag0 - (step1 * 360);

thetag = thetag0 + (360.98564724 * (time / 24));

if thetag > 360
    thetag = thetag - 360;
end

if direction == 'W' || direction == 'w'
    longitude = (longitude * (-1));
end

theta = thetag + longitude;

if theta > 360
    theta = theta - 360;
end

%converting to ECI
Rz = [cosd(theta), sind(theta), 0; -sind(theta), cosd(theta), 0; 0, 0, 1];
Ry = [sind(latitude), 0 ,-cosd(latitude); 0, 1, 0; cosd(latitude), 0, sind(latitude)];

D = Ry * Rz;

inverse = inv(D);

rXYZ = inverse * transpose(rh);
rXYZ = transpose(rXYZ);

omega = [0,0,0.05883];

rhoDot = inverse * transpose(rhoDot);
rhoDot = transpose(rhoDot);

vXYZ = rhoDot + (cross(omega, rXYZ));

y = dot(rXYZ, vXYZ);

r = norm(rXYZ);
v = norm(vXYZ);

fprintf('\nPosition and Velocity Vectors in ECI\n')
fprintf('\nrXYZ = %fx + %fy + %fz DU\n', rXYZ)
fprintf('vXYZ = %fx + %fy  %fz DU\\TU\n', vXYZ)

%specific angular momentum
h = cross(rXYZ, vXYZ);
h1 = norm(h);

z = [0,0,1];
n = cross(z, h);
n1 = norm(n);

%eccentriciy in Canonical Units
eDU = ((v^2 - (1 / r)) * rXYZ) - ((dot(rXYZ,vXYZ)) * vXYZ);
e = norm(eDU);

%Semi-Latus Rectum
p = h1^2;

%Semi-major axis
a = p / (1 - e^2);

%argument of periapsis
w = acosd(dot(n,eDU) / (n1 * e));

if eDU(3) > 0 && w > 180
        w = 360 - w;
else if eDU(3) < 0 && w < 180
                w = 360 - w;
        end
end

%longitude of the ascending node
omega = acosd(n(1) / norm(n));

if n(2) > 0 && omega > 180
        omega = 360 - omega;
else if n(2) < 0 && omega < 180
                omega = 360 - omega;
        end
end

%inclination
i = acosd(h(3) / h1);

%true anomaly
ta = acosd((dot(eDU, rXYZ)) / (e * r));

if dot(rXYZ, vXYZ)
ta = 360 - ta;
end

fprintf('-----------------------------------------------\n\n')


fprintf('Orbital Elements:\n\n')
fprintf('Eccentricity: %.3f\n', e)
fprintf('Semi-Latus Rectum: %.3f DU\n', p)
fprintf('Semi-Major Axis: %.3f DU\n', a)
fprintf('Argument of Periapsis: %.2f degrees\n', w)
fprintf('Longitude of the Ascending Node: %.2f degrees\n', omega)
fprintf('Inclination: %.2f degrees\n', i)
fprintf('True Anomaly: %.2f degrees\n', ta)

a1 = a * 6378;
T = (2 * pi * (a1 ^ (3/2))) / sqrt(mu);
hours = floor(T / 3600);
fprintf('Orbital Period: %.0f seconds\n', T)
minutes = floor((T - (3600 * hours))/ 60);
seconds = (((T - (3600 * hours))/ 60) - minutes) * 60;

fprintf('-----------------------------------------------\n')

fprintf('\nConfirm Object Is In Tundra Orbit\n')
fprintf('\nEccentricity Between 0.2 and 0.3\n')
fprintf('Eccentricity: %.2f\n', e)

if e >= 0.2 && e <= 0.3
    fprintf('Eccentricity Confirmed\n\n')
    check(1) = 1;
else
    fprintf('Eccentricity: NEGATIVE\n\n')
    check(1) = 0;
end

fprintf('Inclination Approximately 63.40 Degrees\n')
fprintf('Inclintation: %.2f\n', i)

if i >= 60 && i <= 65
    fprintf('Inclination Confirmed\n\n')
    check(2) = 1;
else
    fprintf('Inclination: NEGATIVE\n\n')
    check(2) = 0;
end

fprintf('Orbital Period: 1 Sidereal Day\n')
fprintf('Approximate Orbital Period: %.0f hours %.0f minutes %.0f seconds\n', hours, minutes, seconds)

if T > 86064 && T < 86264
    fprintf('Orbital Period Confirmed\n\n')
    check(3) = 1;
else
    fprintf('Orbital Period: NEGATIVE\n\n')
    check(3) = 0;
end


if mean(check) == 1
    fprintf('Target Confirmed: SiriusXM Satellite\n')
else
    fprintf('TARGET NOT CONFIRMED. NOT IN TUNDRA ORBIT\n')
end


E0 = acos((e + cosd(ta)) / (1 + (e * cosd(ta))));
E0D = rad2deg(E0);
if E0D < 180
    E0 = (2 * pi) - E0;
end

M0 = (E0 - (e * sin(E0)));
M0D = rad2deg(M0);
n0 = sqrt((mu) / (a1^3));

dtp = (1 / n0) * M0;

tof = (7 * 3600) + (25 * 60);

k = floor((dtp + tof) / T);

M = (n0 * (tof)) - ( 2 * pi * k) + M0;

Ei = 0.5;
ip = 1;
 
while ip < 10
    
    table(ip,1) = ip;
    table(ip,2) = Ei;
    
    dM = M - (Ei - (e * sin(Ei)));
    dMde = 1 - (e * cos(Ei));
    
    table(ip,3) = dM;
    table(ip,4) = dMde;
    
    Ei = Ei + (dM / dMde);
    table(ip,5) = Ei;
    ip = ip + 1;
end

ta2 = acosd((cos(Ei) - e) / (1 - (e * cos(Ei))));


fprintf('-----------------------------------------------\n')
fprintf('\nThe True Anomaly of the SiriusXM Satellite after %.0f hours and %.0f minutes after loss of signal is %.2f degrees\n', htof, mtof, ta2)

muCU = 1;

r2 = p / (1 + (e * cosd(ta2)));
rw(1,1) = (r2 * cosd(ta2));
rw(2,1) = (r2 * sind(ta2));
rw(3,1) = 0;


vw(1,1) = sqrt(muCU / p) * (-sind(ta2));
vw(2,1) = sqrt(muCU / p) * (e + cosd(ta2));
vw(3,1) = 0;

rW = [cosd(w) -sind(w) 0; sind(w) cosd(w) 0; 0 0 1 ];
rI = [1 0 0; 0 cosd(i) -sind(i) ; 0 sind(i) cosd(i)];
rOmega = [cosd(omega) -sind(omega) 0; sind(omega) cosd(omega) 0; 0 0 1 ];

rotation = rOmega*rI*rW;

rECI = rotation * rw;
vECI = rotation * vw;

fprintf('\nEarth-Centric Equatorial Cartesian Coordinates:\n')
fprintf('\nPostion Vector: %fx  %fy + %fz DU\n',rECI(1,1), rECI(2,1), rECI(3,1))
fprintf('Velocity Vector: %fx + %fy + %fz DU/TU\n\n', vECI(1,1), vECI(2,1), vECI(3,1))



