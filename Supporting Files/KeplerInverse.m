E0 = acos((e + cosd(ta)) / (1 + (e * cosd(ta))));
E0D = rad2deg(E0);

M0 = (E0 - (e * sin(E0)));
M0D = rad2deg(M0);
n0 = sqrt((mu) / (a1^3));

dtp = (1 / n0) * M0;

tof = (htof * 3600) + (mtof * 60);

k = floor((dtp + tof) / T);

M = (n0 * (tof)) - ( 2 * pi * k) + M0;

Ei = 0.5;
ip = 1;
 
while ip < 10
    
    dM = M - (Ei - (e * sin(Ei)));
    dMde = 1 - (e * cos(Ei));
    
    Ei = Ei + (dM / dMde);
    ip = ip + 1;
end

ta2 = acosd((cos(Ei) - e) / (1 - (e * cos(Ei))));

fprintf('-----------------------------------------------\n')
fprintf('\nThe True Anomaly of the SiriusXM Satellite after %.0f hours and %.0f minutes after loss of signal is %.2f degrees\n', htof, mtof, ta2)