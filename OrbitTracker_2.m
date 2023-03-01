%this is the slimmed down version of the OrbitTracker script. 
%It is set up to take different inputs for:
%Range
%Range Rate
%Radar Site Location and Time of Observation
%Time of flight

%all of the scripts here are comprised of the functions in OrbitTracker
%here are the lines that each script corresponds to in OrbitTracker for
%reference
%Sidereal: Lines 40 - 89
%convertToECI: Lines 92 - 115
%orbitalElements: Lines 118 - 178
%tundraCheck: Lines 180 - 221
%KeplerInverse: Lines 224 - 254
%perifocalToECI: Lines 256 - 279

%this scripts inputs the range and range rate vectors (IN CANONICAL UNITS)
%of a satellite in a Tundra orbit around Earth. It also inputs the location
%of the radar site, time zone of the observation, time of the observation, 
%and time of flight between two points.
%With these inputs, the script:

%1. Converts the intial range and range rate vectors to ECI position and
%velocity vectors.

%2. Calculates the orbital elements of the satellite.

%3. Confirms the satellite is in a tundra orbit

%4. Determines the true anomaly at the new point.

%5. Calculates the ECI position and velocity at this new point.

%Anthony D'Amico
%Febuary 25th, 2022

%gravitational parameter of Earth
mu = 3.986004 * 10^5;

%inputs range in canonical units
range  = 'Range: ';
rho = input(range);

%input the range rate in canonical units
rangeRate = 'Range Rate: ';
rhoDot = input(rangeRate);

fprintf('\n')

%convert the given time to Sidereal Time
%can take time from any continental US time zones or UT (in 24 hour format)
run Sidereal

latitude = 'Latitude: ';
latitude = input(latitude);

dir = 'N or S: ';
NoS = input(dir, 's');

fprintf('\nTime of Flight: \n')

prompt = 'Hours: ';
htof = input(prompt);

prompt = 'Minutes: ';
mtof = input(prompt);

fprintf('-----------------------------------------------\n\n')

%convert range and range rate to ECI postion and velocity vectors
run convertToECI

%using r and v to determine orbital elements
run orbitalElements

%confirms that the orbital elements are consistent with a tundra orbir
run tundraCheck

fprintf('-----------------------------------------------\n\n')

%find the new True Anomaly after the given time of flight
run KeplerInverse

%uses orbital elements and new TA to calculate new perifocal r and v
%vectors, then converts them to ECI
run perifocalToECI