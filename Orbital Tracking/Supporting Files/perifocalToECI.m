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