function [Re_G,Re_L]=Reynold(Fvec,T,P,d_P,CSA,voidage)
%Calculates Gas and Liquid Reynolds Numbers
%Inputs:
%       Fvec: Component Molar Flowrate [mol/s]
%       Q_G,Q_L:Gas and liquid volume flowrates [m3/s]
%       d_P:Catalyst Particle Diameter [m]
%       CSA:ReactorCross Sectional Area
%       voidage: bed porosity
%--------------------------------------------------------------------------
[Q_G,Q_L]=VolumeFlow(Fvec,T,P);%returns volume flows
[visc_gas,visc_liq]= viscosity(T,P,Fvec);%returns kinematic viscosity
U_G=Q_G/CSA;%gas superficial velocity [m/s]
U_L=Q_L/CSA;%liquid superficial velocity [m/s]
Re_G=U_G*d_P/(visc_gas*(1-voidage));
Re_L=U_L*d_P/(visc_liq*(1-voidage));
end