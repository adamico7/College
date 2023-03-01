function [distance] = TBCC(FMFR,fuelMASS)

%README
%this code analyizes the preformance of a Turbine Based Combined Cycle
%engine that powers an aircraft. The three distinct operating regimes are
%turbojet, turbojet with afterburner, and ramjet. These three are
%indpendent of each other and only switch when the TSFC from the next phase
%is lower than the previous phase. This allows the TBCC to be the most
%efficient it can be. The turbojet and ramjet are both Brayton cycles.

%The TBCC function takes two inputs. The first is the fuel mass flow rate
%(in kg/s) and the second is the amount of fuel the aircraft has at takeoff
%(in kg). The function returns the total range of the aircraft for the
%given parameters and greats plots detailing the performance of the
%aircraft and individual engines over the flight duration. 

%ASSUMPTIONS:
%Ramjet and Turbojet are both ideal and follow ideal Brayton cycles.
%The nozzles are isentropic and adiabatic.
%The air to fuel ratio is stoichiometric, and since the fuel mass flow rate
%is constant, so is the air mass flow rate.
%Air is an ideal gas.




h = 0;
i = 1;

p = 101.325; %sea level pressure
gamma = 1.285; %mixture
R = 8.31446;
Cp = 1.307; %mixture

%air to fuel ratio
AFRatio = (4.76 * 29 * 15.5) / (142);

massTJ = 1740; %GE-79 mass
massRJ = 459;  %XRJ47 mass
mass(i) = fuelMASS + massTJ + massRJ; %total aircraft mass at takeoff

timeOfFlight = 0;
throwAway = 0;

check = 0;
checks = 0;
peach = 0;

Range(i) = 0;

while fuelMASS > 0
    
    H(i) = h;
    
    %atmospheric model generator
    
   if h < 11.4
       %troposphere
       
       a = (-56.9 - 15) / 11; %temperature lapse rate
       T(i) = 288.15 + (a * h); %temperature in troposphere
       p(i) = 101.325 * (T(i) / 288.15)^(5.256); %pressure in troposphere
   
   else
       %tropopause
       T(i) = 213.55; %constant temperature in tropopause
       p(i) = 22.56 * exp(1).^(1.73 - 0.000157*(h*1000)); %tropopause pressure
   
   end
   
       rho(i) = p(i) / (R * T(i)); %density based on ideal gas law
       sos(i) = sqrt(gamma * R * T(i) / 0.029); %speed of sound 
   
   
       ma(i) = FMFR * AFRatio; %air mass flow rate
   
       f = FMFR / ma(i); %fuel to air ratio
   
   %for take off
    if h <= 0.1 && throwAway == 0
    
       
    mass(i) = mass(i);
    T3(i) = 1210;
    T4(i) = 1210;
    T0(i) = T(i);
    v1(i) = 0; %not moving for takeoff
    range = 0;
    
    
    v4(i) = 1345.03819363999; %exhaust velocity for 1210 K
    
    SFC(i) = ((1 + f) * v4(i)) - v1(i); 
    thrust(i) = SFC(i) * ma(i); 
    
    acceleration = thrust(i) / mass(i);
    v1(i) = v1(i) + acceleration; 
    
    TSFC(i) = f / SFC(i);
    
   
    end
    
    TSFCcheck = 0;

    while TSFCcheck < 1 &&  i > 1

    %subsonic cruise 
        if MACH(i - 1) < 1
    
            mass1(i) = mass(i-1) - (FMFR);
    
            T31(i) = 1210;
            T41(i) = T31(i) / ((pressureRatio(i-1)) ^ ((gamma - 1) / gamma));
            T041(i) = T41(i) * (1 + (((gamma - 1) / 2) * MACH(i-1) ^ 2));
    
            v41(i) = 1345.03819363999;
    
            blah = v1(i - 1);
    
            SFC1(i) = ((1 + f) * v41(i)) - blah;
    
    
            if v1(i-1) == v41(i)
                thrust1(i) = 0;
            else
                thrust1(i) = SFC1(i) * ma(i);
            end
    
            acceleration = thrust1(i) / mass1(i);
    
            if (acceleration + v1(i - 1) > v41(i))
                v11(i) = v41(i);
            else
                v11(i) = v1(i - 1) + acceleration;
            end
    
            if SFC1(i) == 0
                TSFC(i) = 0;
            else
                TSFC1(i) = f / SFC1(i);
            end
    
    
            range = 0;
    
            if h > 20
                range1 = v11(i) * ((mass1(i) * 9.81) / (ma(i) * v1(i-1))) * log(mass(i-1) / mass1(i)) * (ma(i) * v1(i-1) / FMFR) / 9.81;
                range1 = range1 / 1000;
            else
                range1 = 0;
            end
    
            timeOfCruise = i;
            M = mass1(i);
            pepe = i;
            speed = v11(i);
            press = pressureRatio(i-1);
    
        else %afterburner preformance
            p62 = p2(i-1);
            p52 = p62;
            T6 = 1210;
            T7(i) = T6 - (0.5 * (v4(i-1)^2) / (Cp * 1000));
            T07(i) = T7(i) * (1 + (((gamma - 1) / 2) * MACH(i-1) ^ 2));
            p62 = p0(i - 1) *((T6 / T7(i))^(gamma/(gamma - 1)));
        
            mass1(i) = mass(i-1) - (FMFR);
        
            T31(i) = 1210;
            T41(i) = T31(i) / ((pressureRatio(i-1)) ^ ((gamma - 1) / gamma));
            T041(i) = T41(i) * (1 + (((gamma - 1) / 2) * MACH(i-1) ^ 2));
        
        
            v41(i) = v1(i - 1) * sqrt((T7(i-1)) / T0(i-1));
            if v41(i) == 0
                v41(i) = v41(i-1);
            end
        
       
            SFC1(i) = ((1 + f) * v41(i-1)) - v1(i-1);
        
            if SFC1(i) < 0
                SFC1(i) = SFC1(i - 1);
            elseif SFC1(i) < (SFC1(i-1) / 10)
                SFC1(i) = SFC1(i-1);
            end
        
            if v41(i) < (v41(i-1) / 3)
                v41(i) = v41(i-1);
            end
        
            thrust1(i) = SFC1(i) * ma(i);
    
            acceleration = thrust1(i) / mass1(i);
            v11(i) = v1(i - 1) + acceleration;
        
            TSFC1(i) = f / SFC1(i);
        
            if h >= 20
                range1 = v11(i) * ((mass1(i) * 9.81) / (ma(i) * v1(i-1))) * log(mass(i-1) / mass1(i)) * (ma(i) * v1(i-1) / FMFR) / 9.81;
                range1 = range1 / 1000;
            
            else
                range1 = 0;
            end
            
            
        end
    

        
        %ramjet
        T33(i) = 900;
        T43(i) = T33(i) / ((pressureRatio(i-1)) ^ ((gamma - 1) / gamma));
        T043(i) = T43(i) * (1 + (((gamma - 1) / 2) * (MACH(i-1) ^ 2)));
        
        mass3(i) = mass(i-1) - (FMFR);
     
        v43(i) = v1(i - 1) * sqrt(T043(i-1) / T0(i - 1));
        
        SFC3(i) = ((1 + f) * v43(i-1)) - v1(i-1);
        
        thrust3(i) = SFC3(i) * ma(i);
        
        acceleration(i) = thrust3(i) / mass3(i);
        
        if v1(i - 1) + acceleration(i) > v43(i)
            v13(i) = v43(i);
        else
            v13(i) = v1(i - 1) + acceleration(i);
        end
        
        if SFC3(i) == 0
            TSFC3(i) = 0;
        else
            TSFC3(i) = f / SFC3(i);
        end
        
        
        if h >= 20
        
            range3 = v13(i) * ((mass3(i) * 9.81) / (ma(i) * v1(i-1))) * log(mass(i-1) / mass3(i)) * (ma(i) * v1(i-1) / FMFR) / 9.81;
            range3 = range3 / 1000;
        
        else
            range3 = 0;
        end

        TSFCcheck = TSFCcheck + 1;
        
        if TSFC1(i) > TSFC3(i) && peach < 1 && i > 10
            error = i;
            peach = 1;
        end
        
    end
    
    if i > 1
    
        if  TSFC1(i) > TSFC3(i) && TSFC3(i) > 0
            TSFC(i) = TSFC3(i);
            mass(i) = mass3(i);
            T3(i) = T33(i);
            T4(i) = T43(i);
            T04(i) = T043(i);
            v4(i) = v43(i);
            SFC(i) = SFC3(i);
            thrust(i) = thrust3(i);
            v1(i) = v13(i);
            range = range3;
    
            check = check + 0.5;
            check2 = 2;
    
    
        elseif check ~= 2 
            
            TSFC(i) = TSFC1(i);
            mass(i) = mass1(i);
            T3(i) = T31(i);
            T4(i) = T41(i);
            T04(i) = T041(i);
            v4(i) = v41(i);
            SFC(i) = SFC1(i);
            thrust(i) = thrust1(i);
            v1(i) = v11(i);
            range = range1;
    
            long = i;
            
        else
            TSFC(i) = TSFC3(i);
            mass(i) = mass3(i);
            T3(i) = T33(i);
            T4(i) = T43(i);
            T04(i) = T043(i);
            v4(i) = v43(i);
            SFC(i) = SFC3(i);
            thrust(i) = thrust3(i);
            v1(i) = v13(i);
            range = range3;

        end
        
    end
 
    MACH(i) = v1(i) / sos(i);
    MACH2(i) = v4(i) / sos(i);

    T0(i) = T(i) * (1 + (((gamma - 1) / 2) * MACH(i) ^ 2));
    p0(i) = p(i) * ((T0(i) / T(i)) ^ (gamma / (gamma - 1)));
    T2(i) = T0(i);
    p2(i) = p0(i);
    p3(i) = p2(i); %constant pressure heat addition
    p4(i) = (p3(i) * p(i)) / p0(i);

    pressureRatio(i) = p2(i) / p(i);

    Wexp(i) = Cp * (T3(i) - T4(i));
    Wcomp(i) = Cp * (T2(i) - T(i));
    thermEff(i) = 1 - (1 / ((p2(i) / p(i))^((gamma - 1) / gamma)));
    propEff(i) = 2 / (1 + (v4(i) / v1(i)));
    overallEff(i) = thermEff(i) * propEff(i);

    propEff(i) = propEff(i) * 100;
        if propEff(i) > 100
            propEff(i) = propEff(i - 1);
        end
    thermEff(i) = thermEff(i) * 100;
    overallEff(i) = overallEff(i) * 100;
    
    if i > 2
        if overallEff(i) > (overallEff(i-1) + 10)
            overallEff(i) = overallEff(i - 1);
    
        end
    end

    fuelMASS = mass(i) - massTJ - massRJ;

    if h < 20
        h = h + (v1(i) / 1000);
    
    end
    
    i = i + 1;
   
    timeOfFlight(i) = timeOfFlight(i-1) + 1;
    Range(i) = range + Range(i - 1);
   
    throwAway = throwAway + 1;

end

%error section smoothes graphs as the engines transition
TSFC(error) = TSFC(error - 2);
TSFC(error + 1) = TSFC(error);
TSFC1(error) = TSFC(error - 1);
TSFC1(error + 1) = TSFC(error);

thermEff(error) = thermEff(error - 2);
thermEff(error + 1) = thermEff(error);

propEff(error) = propEff(error - 2);
propEff(error - 1) = propEff(error);
propEff(error + 1) = propEff(error);

overallEff(error) = overallEff(error - 2);
overallEff(error - 1) = overallEff(error);
overallEff(error + 1) = overallEff(error);

figure;
plot(timeOfFlight(1:length(thrust)), thrust)
hold on
ylabel('Thrust')
yyaxis right
plot(timeOfFlight(1:length(thrust)), v4)
plot(timeOfFlight(1:length(thrust)), v1)
yline((4 * sos(i-1)), 'r', 'Mach 4');
yline(1 * sos(i-1), 'b', 'Mach 1');
hold off
xlabel('Time of Flight (s)')
ylabel('Velocity (m/s)')
legend('Thrust','Exit Velocity', 'Inlet Velocity')

figure;
plot(Range(1:length(thrust)), thrust)
hold on
ylabel('Thrust')
yyaxis right
plot(Range(1:length(thrust)), v4)
plot(Range(1:length(thrust)), v1)
yline((4 * sos(i-1)), 'r', 'Mach 4');
yline(sos(i-1), 'b', 'Mach 1');
hold off
xlabel('Range (km)')
ylabel('Velocity (m/s)')
legend('Thrust','Exit Velocity', 'Inlet Velocity')

figure;
plot(timeOfFlight(1:length(thrust)), overallEff)
hold on
plot(timeOfFlight(1:length(thrust)), thermEff)
plot(timeOfFlight(1:length(thrust)), propEff)
hold off
xlabel('Time of Flight (s)')
ylabel('Efficiency (%)')

legend('Overall Efficiency','Thermal Efficinecy','Propulsion Efficinecy')

figure;
plot(timeOfFlight(1:length(thrust)), TSFC)
xlabel('Time of Flight (s)')
ylabel('TSFC')

distance = Range(length(Range));

figure;
plot(timeOfFlight(1:length(thrust)), H)
xlabel('Time of Flight (s)')
ylabel('Altitude (km)')

figure;
plot(timeOfFlight(107:750), TSFC1(107:750))
hold on
plot(timeOfFlight(107:750), TSFC3(107:750))
hold off
xlabel('Time of Flight (s)')
ylabel('TSFC')
legend('TSFC (Turbojet)','TSFC (Ramjet)')

distance = Range(length(Range));


end

