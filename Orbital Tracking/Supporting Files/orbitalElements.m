r = norm(rXYZ);
v = norm(vXYZ);

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