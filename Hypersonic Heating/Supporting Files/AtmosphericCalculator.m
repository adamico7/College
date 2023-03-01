%% Documentation:

%This function computes the atmospheric properties on Earth for any
%altitude below 150km. Input geometric altitude in km above sea level.
%Outputs temperature in Kelvin, pressure in Pa, density in kg/m^3, and
%enthalpy in J/kg-K.

%Validated against 1976 standard atmosphere.

function [T1,p1,rho1,a1] = AtmosphericCalculator(altitude)
%altitude = altitude*1000; %convert to km
mw = 28.85;

if altitude <= 86
    h = (altitude .* 6371) / (altitude + 6371);
    
    if h <= 11
        lhi = -6.5;
        tmi = 288.15;
        pi = 101325;
        hi = 0;
    elseif h <= 20
        lhi = 0;
        tmi = 216.65;
        pi = 22631.95;
        hi = 11;
    elseif h <= 32
        lhi = 1;
        tmi = 216.65;
        pi = 5474.79;
        hi = 20;
    elseif h <= 47
        lhi = 2.8;
        tmi = 228.65;
        pi = 868.01;
        hi = 32;
    elseif h <= 51
        lhi = 0;
        tmi = 270.65;
        pi = 110.9;
        hi = 47;
    elseif h <= 71
        lhi = -2.8;
        tmi = 270.65;
        pi = 66.94;
        hi = 51;
    elseif h <= 84.855
        lhi = -2;
        tmi = 214.65;
        pi = 3.956;
        hi = 71;
    end
    if lhi ~= 0
        T1 = tmi + lhi * (h-hi);
        p1 = pi *(tmi/T1) .^ ((9.806 * 28.9664)/(8.31432 * lhi));
        rho1 = (p1 * 28.9664)/(8314*T1);
    elseif lhi == 0
        T1 = tmi + lhi * (h-hi);
        p1 = pi * exp((-9.806*28.9664*(h-hi))/(8.314*T1));
        rho1 = (p1 * 28.9664)/(8314*T1);
    end
else
    if altitude < 100
        lzi = 1.6481;
        tmi = 186.945;
        pi = 0.344180;
        zi = 86;
        mw = 28.9644;
    elseif altitude < 110
        lzi = 5;
        tmi = 210.65;
        pi = 0.029073;
        zi = 100;
        mw = 28.88;
    elseif altitude < 120
        lzi = 10;
        tmi = 260.65;
        pi = 0.006801;
        zi = 110;
        mw = 28.56;
    elseif altitude < 150
        lzi = 20;
        tmi = 360.65;
        pi = 0.002247;
        zi = 120;
        mw = 28.08;
    end
    R = 8314/mw;
    b = 3.31*10^-7;
    powP = -((9.806/(R*(lzi/1000)))*(1+b*((tmi/(lzi/1000)) - zi*1000)));
    ti = (tmi * mw)/28.9664;
    rhoI = ((pi * mw)/(8314*ti));
    powRho = -((9.806/(R*(lzi/1000)))*(((R*(lzi/1000))/9.806) + 1 + b* ((tmi/(lzi/1000)) - zi*1000)));
    term1 = ((tmi + lzi*(altitude - zi))/tmi);
    term2 = exp(((9.806*b)/(R*(lzi/1000))*1000*(altitude-zi)));
    p1 = pi * term1 .^ powP * term2;
    rho1 = rhoI * term1 .^ powRho .* term2;
    T1 = p1 * mw/ (8314 * rho1);
end
        
a1 = sqrt(1.4*287*T1);
cpAir = 1004.5;
h1 = T1 * cpAir;
end