function dP=Holub(Fvec,T,P,d_P,CSA,voidage,cat_density,dzdW)
%Calculates Pressure Drop Across Reactor
%--------------------------------------------------------------------------
%External Functions and Constants
E1=150;E2=1.75;%Ergun Parameters-Perfect Spheres assumed
[Re_G,Re_L]=Reynold(Fvec,T,P,d_P,CSA,voidage);
[Ga_G,Ga_L]=Galileo(T,P,Fvec,voidage,d_P);
[~,liquid_holdup]= Rao(d_P,cat_density,Fvec,T,P,voidage);
[~,mass_density] = density(T,P);
g=9.81;%gravitational acceleration [m/s2]
%-------------------------------------------------------------------------- 
%Density Calculations
%Gas Density and Mass Flow
rho_G=mass_density(1);%Hydrogen Density

%Liquid Density and Mass Flow
Fvec_L=Fvec(2:6);%liquid molar flows [kmol/s]
RFM=[92.09382;74.08;76.09;18.01528;60.09502];%[G;A;P;POH][kg/kmol]
Fvec_mass=bsxfun(@times,Fvec_L,RFM);%mass flow vector [kg/s]
m_L=sum(Fvec_mass);%liquid mass flow [kg/s]
x_i=Fvec_mass/m_L;%mass fraction vector
rho_i=mass_density(2:6);%liquid mass density vector [kg/m3]
v_molar=1./rho_i;%liquid specific volume vector [m3/kg]
v_L=dot(x_i,v_molar);%liquid mixture specific volume [m3/kg]
rho_L=1/v_L;% Liquid mixture density [kg/m3]
%--------------------------------------------------------------------------
%Dimensionless Pressure force calculation
%Gas
a_G=(voidage/(voidage-liquid_holdup))^3;
b_G=E1*Re_G/Ga_G;
c_G=E2*(Re_G^2)/Ga_G;
PSI_G=a_G*(b_G+c_G);
%--------------------------------------------------------------------------
%Dimensionless Pressure force calculation
%Liquid
a_L=(voidage/liquid_holdup)^3;
b_L=E1*Re_L/Ga_L;
c_L=E2*(Re_L^2)/Ga_L;
PSI_L=a_L*(b_L+c_L);
%--------------------------------------------------------------------------
dP_Gdz=rho_G*g*(PSI_G-1);%gas ?P wrt reactor length
dP_G=dP_Gdz*dzdW;%gas ?P wrt catalyst weight
dP_Ldz=rho_L*g*(PSI_L-1);%liquid ?P wrt reactor length
dP_L=dP_Ldz*dzdW;%liquid ?P wrt catalyst weight
%dP vector
dP=[dP_G;dP_L];
end

