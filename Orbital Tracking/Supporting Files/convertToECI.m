Rz = [cosd(theta), sind(theta), 0; -sind(theta), cosd(theta), 0; 0, 0, 1];
Ry = [sind(latitude), 0 ,-cosd(latitude); 0, 1, 0; cosd(latitude), 0, sind(latitude)];

rh = rho + [0,0,1];

D = Ry * Rz;

inverse = inv(D);

rXYZ = inverse * transpose(rh);
rXYZ = transpose(rXYZ);

omega = [0,0,0.05883];

rhoDot = inverse * transpose(rhoDot);
rhoDot = transpose(rhoDot);

vXYZ = rhoDot + (cross(omega, rXYZ));


r = norm(rXYZ);
v = norm(vXYZ);

fprintf('\nPosition and Velocity Vectors in ECI\n')
fprintf('\nrXYZ = %fx + %fy + %fz DU\n', rXYZ)
fprintf('vXYZ = %fx + %fy  %fz DU\\TU\n', vXYZ)