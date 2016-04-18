function [liquid_sat,liquid_holdup]= Rao(d_P,cat_density,Fvec,T,P,voidage)
%Rao: Calculates liquid saturation=volume of liquid per unit void volum
%Neglects surface tension effects
%Particles assumed large enough for this assumption
%Liquid Holdup:
%Volume of liquid per bed volume is subsequently calculated from liquid
%saturation: (liquid holdup)=voidage*(liquid saturation)
%Inputs:
%mass_density vector:[rho_H;rho_G;rho_A;rho_P;rho_W;rho_POH]
%                    calculated using density(T[K],P[Pa]) function
%d_P:particle diameter [m]
%Fvec:component molar flow vector [mol/s]
Fvec=Fvec/1000;%[kmol/s]
%cat_density:catalyst density [g/m3]
cat_density=cat_density/1000;%[kg/m3]
%voidage: catalyst voidage
%-------------------------------------------------------------------------- 
%Density Calculations
%Gas Density and Mass Flow
[~,mass_density] = density(T,P);
rho_H=mass_density(1);%Hydrogen Density
RFM_H=2;%kg/kmol
m_H=rho_H*Fvec(1)*RFM_H;%mass density*molar flowrate*RFM

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
%specific surface area
a_s=6/(cat_density*d_P);
%Liquid Saturation
c=0.4;
b=0.27;%assuming pulse flow
X=(m_L/m_H)*sqrt(rho_H/rho_L);%lockhartmartinelli
liquid_sat=c*X^b*a_s^(1/3);
liquid_holdup=liquid_sat*voidage;
end

