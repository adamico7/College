function [Cl,Cd, q_x, qs] = hypersonics(x,y, M1, h, aoa)

run bodyGeometry.m

format long g

[T1,P1,rho1,a] = AtmosphericCalculator(h);
v = M1 * a;

rho2 = (2.4 * rho1 * M1^2)/(2 + (0.4 * M1^2));
p_ratio = 1 + (2.8 * (M1^2 - 1)/2.4);
p2 = P1 * p_ratio;

T21 = T1 * p2 * rho1 / rho2 / P1;

M2 = sqrt((1 + (0.2 * M1^2))/((1.4*M1^2)-0.2));

a2 = sqrt(1.4 * 287 * T21);
%u2 = M2 * a2;

h1 = 1004.5 * T1;
h1 = h1/1000;

tol = 1 * 10^-6;
ep = rho1/rho2;
erreps = 1;
i = 1;
ep(i) = 0.1;

 while (erreps >= tol)
    
    p2i = P1 + rho1 * v^2 * (1 - ep(i));
    h2i = h1 + (v^2 / 2) * (1 - ep(i)^2);
    h2i = h2i / 1000;

    p2i = p2i / 101325;

    T2j = clib.EquilFlow.EquilFlow.Composition(2,p2i,h2i);
    
    p2i = p2i * 101325;

    
    dim = length(T2j);
    
    MWmix = T2j(dim);

    R = 8314 / MWmix;
    
   
    rho2i = p2i/(R * T2j(10));

    ep(i + 1) = rho1 / rho2i;
    erreps = abs((ep(i)-ep(i+1))/ep(i));

    i = i + 1;

 end

 u2 = ep(i) * v;


Cpmax = 2 - ep(i);
long = length(x);

i = 1;

T2 = T2j(10);
%T2 = T21;

while i < long

    delta(i) = atand((y(i + 1) - y(i))/(x(i + 1) - x(i)));
    theta(i) = 90 - delta(i);

    rn = 0.03048;

    if i > 1

        CAi = (1 - (rn/y(long))^2 * cosd(delta(i))^2) * ((Cpmax * sind(delta(i))^2 * cosd(aoa)^2) +  (Cpmax/2 * cosd(delta(i))^2 * sind(aoa)^2));
        CA(i) = CAi * ((y(i+1)^2 - y(i)^2) / y(long)^2);
        
        
        CNi = Cpmax * (1 - ((rn/y(long))^2 * cosd(delta(i))^2)) * (cosd(delta(i))^2 * sind(aoa) * cosd(aoa));
        CN(i) = CNi * ((y(i+1)^2 - y(i)^2) / y(long)^2);

    else

        delta(2) = atand((y(3) - y(2))/(x(3) - x(2)));
        CA(i) =  (rn/y(long))^2 * (Cpmax/2 * (1 - sind(delta(2))^4));
        CN(i) = 0;


    end

    i = i + 1;

end

theta(1) = theta(2);

CA_value = sum(CA);
CN_value = sum(CN);

Cl = CN_value * cosd(aoa) - CA_value * sind(aoa);
Cd = CN_value * sind(aoa) + CA_value * cosd(aoa);

if M1 == 10 && h == 18.28

    T2j = clib.EquilFlow.EquilFlow.Composition(2,p2i,h2i);
    
    
    T = T2;
    Cp = (8.31935 * 10^-12 * T^3) - (8.62469 * 10^-8 * T^2) + (3.14166 * 10^-4 * T) + 0.901439;
    Cp = Cp * 1000;
    [k, mu] = complicated(T2, p2i / 101325);
    
   
    rho2 = rho2i;
    Pr = mu * Cp / k;
    

    p_inf = P1;
    P_wall = p2i;
    R = 8314 / T2j(14);
    T_wall = 300;
    Cp_air = 1004.5;
    h01 = Cp * T2;
    h0 = T2j(11) * 1000;
    Tw = T_wall;

    hw = Cp_air * T_wall;
    i = 1;

    while i < long
    
    if x(i) == 0
        
        dUe = (1 / rn) * sqrt((2*(P_wall - p_inf))/rho2i);
        rho_wall = P_wall / (287 * T_wall);
        [k_wall, mu_wall] = complicated(T_wall, P_wall / 101325);
        hw = Cp_air * T_wall;

        

        qs = 0.76 * Pr^(-0.6) * (rho2 * mu)^(0.4) * (rho_wall * mu_wall)^(0.1) * sqrt(dUe) * (h0 - hw) / 10000;

    else

        T = T2;
        gamma_e = (2.39683 * 10^-19 * T^5) - (3.0436 * 10^-15 *T^4) +(8.89216 * 10^-12 * T^3) + (2.77835 * 10^-8 * T^2) - (1.46939 * 10^-4 * T) + (1.4517);
        a_e = sqrt(gamma_e * (p2i / rho2));
        M_e = u2 / a_e;

        T_ref = (T2) * (1 + 0.032 * M_e^2 + 0.58 * ((Tw / T2) - 1));
        



        [k_star, mu_star] = complicated(T_ref, p2i / 101325);

        T = T_ref;

        Cp_e = (8.31935 * 10^-12 * T^3) - (8.62469 * 10^-8 * T^2) + (3.14166 * 10^-4 * T) + 0.901439;
        Pr_ref = mu_star * Cp_e / k_star * 1000;
        [k2, mu2] = complicated(T2, p2i / 101325);
        Pr = mu2 * Cp / k2;

        i = 1;

        while i < long
            
            if i == 1
                S(i) = rn * deg2rad(theta(i));

            else

                S(i) = sqrt((x(i + 1) - x(i))^2 + (y(i + 1) - y(i))^2);
                S(i) = S(i) + S(i - 1);

            end

            rho_ref = p2i / (R * T_ref);
            Re_ref(i) = rho_ref * u2 * (S(i)) / mu_star;

            if Re_ref(i) <  10^5
                mf = sqrt(3);

                Cf(i) = 0.664 * mf / sqrt(Re_ref(i));
                Ch = 0.332 * (Pr_ref^(-2/3)) * mf / sqrt(Re_ref(i));
                r = sqrt(Pr);


            else
                mf = sqrt(2);

                Cf(i) = 0.0592 * mf / (Re_ref(i))^0.2;
                Ch(i) = 0.5 * Cf(i) * (Pr_ref)^(-2/3);
                r = Pr^(1/3);

                round(Cf, 6);

            end

            Hr = (h2i * 1000) + r * u2^2 / 2;
            tree(i) = rho2 * u2 * Ch(i);
            qx(i) = rho2 * u2 * Ch(i) * (Hr - hw) / 10000;

            i = i + 1;

        end

    end
    i = i + 1;
    end
    q_x = qx;
else
    q_x = 0;

end
end





