function Z=solvePR(T,P)
%--------------------------------------------------------------------------
%Gas Molar Density: For Hydrogen
%Uses Peng-Robinson Equation of State in Cubic form
%LINK: http://chempaths.chemeddl.org/services/chempaths/?q=book/General%20Chemistry%20Textbook/Solids%2C%20Liquids%20and%20Solutions/1409/critical-temperature-and-pressure
R=8.314459848;%J K-1 mol-1
P_C=1.30e6;%Hydrogen Critical Pressure [Pa]
T_C=33.2;%Hydrogen Critical Temperature [K]
T_r=T/T_C;%Hydrogen Reduced Temperature
w=-0.220;%Hydrogen Acentric Factor
a=0.457235*(R^2)*(T_C^2)/P_C;
b=0.077796*R*(T_C^2)/P_C;
k=0.37464+1.54226*w-0.26992*(w^2);
alpha=(1+k*(1-T_r^0.5))^2;
A=a*alpha*P/((R^2)*(T^2));
B=b*P/(R*T);
    function PR=PR(A,B,Z)
        PR=Z^3-(1-B)*Z^2+(A-2*B-3*B^2)*Z-(A*B-B^2-B^3);
    end
    function dPR=dPR(A,B,Z)
        dPR=3*Z^2-2*(1-B)*Z+(A-2*B-3*B^2);
    end
nmax=1e6;
n=1;
tol=1e-6;
difference=tol+1;
if P<4.3e6
    Zold=1;%Initial guess of compressibility factor;
else
    Zold=2;%at P>4.3MPa Zold=1 Newton-Raphson iterates towards wrong root
end

while difference>tol
    if n>nmax
        error('Cannot Calculate Z to specified tolerance')
    else
    n=n+1;
    Znew=Zold-PR(A,B,Zold)/dPR(A,B,Zold);
    difference=abs(Znew-Zold);
    Zold=Znew;
    end
end
Z=Zold;
%--------------------------------------------------------------------------
end