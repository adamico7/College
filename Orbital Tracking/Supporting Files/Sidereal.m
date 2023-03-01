fprintf('Date: \n')
prompt = 'Month: ';
m = input(prompt);

prompt = 'Day: ';
d = input(prompt);

prompt = 'Year: ';
y = input(prompt);

fprintf('\nTime: \n')
prompt = 'Hour: ';
hour = input(prompt);

prompt = 'Minutes: ';
min = input(prompt);

prompt = 'Time Zone: ';
tz = input(prompt, 's');

if tz == 'EST'
    hour = hour + 5;
elseif tz == 'CST'
    hour = hour + 6;
elseif tz == 'MST'
    hour = hour + 7;
elseif tz == 'PST'
    hour = hour + 8;
end
    
time = hour + (min / 60);

fprintf('\nLocation: \n')
prompt = 'Longitude: ';
longitude = input(prompt);

prompt = 'E or W: ';
direction = input(prompt, 's');


J0 = ((367 * y) - (floor(7 * (y + (floor((m + 9)/12))))/4) + floor((275 * m) / 9) + d + 1721013.5);
T0 = (J0 - 2451545) / 36525;



thetag0 = 100.4606184 + (36000.77004 * (T0)) + (0.000387933 * (T0 ^ 2)) + - ((2.583 * 10 ^ -8) * (T0 ^ 3));

thetag = thetag0 + (360.98564724 * (time / 24));
step = floor(thetag / 360);
thetag = thetag - (360 * step);

if direction == 'W' || direction == 'w'
    longitude = (longitude * (-1));
else
    longitude = longitude;
end

theta = thetag + longitude;

if theta > 360
    theta = theta - 360;
end

%fprintf('\nTheta(g) = %f degrees\n', thetag)
%fprintf('Theta = %f degrees\n', theta)