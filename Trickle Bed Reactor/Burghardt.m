function wetting_efficiency=Burghardt(Fvec,T,P,d_P,CSA,voidage)
%Calculates Wetting Efficiency: Fraction of catalyst surface covered by
%liquid phase
%Uses Burghardt et al (1995) correlation
%Inputs:
%       Gas and Liquid Reynold Numbers
%       Liquid Viscosity
%--------------------------------------------------------------------------
[Re_G,Re_L]=Reynold(Fvec,T,P,d_P,CSA,voidage);
[~,visc_liq]= viscosity(T,P,Fvec);
g=9.81;
wetting_efficiency=3.38*(Re_L^0.22)*(Re_G^-0.83)*(d_P*sqrt(g)/visc_liq)^-0.512;
end