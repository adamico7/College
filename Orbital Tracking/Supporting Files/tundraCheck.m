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

fprintf('\n')