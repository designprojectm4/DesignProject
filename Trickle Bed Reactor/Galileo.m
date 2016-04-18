function[Ga_G,Ga_L]=Galileo(T,P,Fvec,voidage,d_P)
%Calculates Galileo Numbers of gas and liquid phases
%Inputs:
%Gas and Liquid Kinematic Viscosities
%[visc_gas,visc_liq]=viscosity(T[K],P[Pa],Fvec[mol/s]) [m2/s]
%d_P: Catalyst Particle diameter [m]
%voidage: bed porosity
%--------------------------------------------------------------------------
[visc_gas,visc_liq]=viscosity(T,P,Fvec);
g=9.81;%gravitational acceleration [m/s2]
Ga_G=g*(d_P^3)*(voidage^3)/((visc_gas^2)*((1-voidage)^2));
Ga_L=g*(d_P^3)*(voidage^3)/((visc_liq^2)*((1-voidage)^2));
end