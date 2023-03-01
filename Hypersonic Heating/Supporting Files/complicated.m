function [mixK, mixViscosity] = complicated(T, p)
format long g

X=clib.EquilFlow.EquilFlow.Composition(1,T,p);


mw(1) = 28.014; %N2
mw(2) = 31.998; %O2
mw(3) = 14.007; %N
mw(4) = 15.999; %O 
mw(5) = 30.006; %NO

X(1) = 0.74423;
X(2) = 0.15558;
X(3) = 2.606 * 10^-5;
X(4) = 0.04123;
X(5) = 0.05887;

Cp(1) = 37.157 / mw(1) * 1000;
Cp(2) = 39.848 / mw(2) * 1000;
Cp(3) = 24.461 / mw(3) * 1000;
Cp(4) = 22.044 / mw(4) * 1000;
Cp(5) = 37.247 / mw(5) * 1000;

Cp(1) = 35.727 / mw(1) * 1000;
Cp(2) = 39.114 / mw(2) * 1000;
Cp(3) = 21.135 / mw(3) * 1000;
Cp(4) = 21.027 / mw(4) * 1000;
Cp(5) = 36.963 / mw(5) * 1000;

sigma(1) = 3.798;
sigma(2) = 3.467;
sigma(3) = 3.298;
sigma(4) = 3.050;
sigma(5) = 3.492;

ek(1) = 71.4;
ek(2) = 106.7;
ek(3) = 71.4;
ek(4) = 106.7;
ek(5) = 116.7;

k1 = 1.380649 * 10^ -23;

Tstar(1) = (1 / ek(1)) * T;
Tstar(2) = (1 / ek(2)) * T;
Tstar(3) = (1 / ek(3)) * T;
Tstar(4) = (1 / ek(4)) * T;
Tstar(5) = (1 / ek(5)) * T;

omega(1) = (1.16145 / (Tstar(1)^0.14874)) + (0.52487 / exp(0.7732 * Tstar(1))) + (2.16178 / (exp(2.43787 * Tstar(1))));
omega(2) = (1.16145 / (Tstar(2)^0.14874)) + (0.52487 / exp(0.7732 * Tstar(2))) + (2.16178 / (exp(2.43787 * Tstar(2))));
omega(3) = (1.16145 / (Tstar(3)^0.14874)) + (0.52487 / exp(0.7732 * Tstar(3))) + (2.16178 / (exp(2.43787 * Tstar(3))));
omega(4) = (1.16145 / (Tstar(4)^0.14874)) + (0.52487 / exp(0.7732 * Tstar(4))) + (2.16178 / (exp(2.43787 * Tstar(4))));
omega(5) = (1.16145 / (Tstar(5)^0.14874)) + (0.52487 / exp(0.7732 * Tstar(5))) + (2.16178 / (exp(2.43787 * Tstar(5))));

i = 1;
while i <= 5

    viscosity(i) = 2.6693 * 10^-5 * ((sqrt(mw(i) * T) / (sigma(i)^2 * omega(i))));
    viscosity(i) = viscosity(i) / 10;
    
    if i == 3 || i == 4
    
        k(i) = 1.9891 * 10^-4 * ((sqrt(T / mw(i)) / (sigma(i)^2 * omega(i))));
        k(i) = k(i) * 418.6;

    end

    if i == 1 || i == 2 || i == 5
        k(i) = viscosity(i) * (Cp(i) + ((5/4) * (8314 / mw(i))));

    end

    i = i + 1;

end

i = 1;

while i <= 5

    j = 1;

    while j <= 5

        psi(i, j) = (1 / sqrt(8)) * (1/(1 + (mw(i) / mw(j)))^(0.5)) * (1 + (viscosity(i)/viscosity(j))^(0.5) * (mw(j) / mw(i))^(0.25))^2;

        j = j + 1;

    end

    i = i + 1;

end

i = 1;

while i <= 5

    j = 1;
    num = (X(i) * viscosity(i));
    num1 = (X(i) * k(i));

    while j <= 5

        den(j) = X(j) * psi(i, j);

        j = j + 1;

    end

    mixtureViscosity(i) = num / sum(den);
    mixtureK(i) = num1 / sum(den);

    i = i + 1;

end

mixViscosity = sum(mixtureViscosity);
mixK = sum(mixtureK);
end