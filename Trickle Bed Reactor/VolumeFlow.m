function[Q_G,Q_L]=VolumeFlow(Fvec,T,P)
%Finds volumetric flows based on molar flows
%   Inputs:
%       Temperature [K]
%       Pressure [Pa]
%       molar density[rho_H;rho_G;rho_A;rho_P;rho_W;rho_POH][kmol/m3]
%       Fvec=molar flow vector [mol/s]
%--------------------------------------------------------------------------
[molar_density,~]=density(T,P);
molar_density=molar_density*1e3;%convert to mol/m3
molar_volume=1./molar_density;%m3/mol
%Gas volume flow(m3/mol)*(mol/s)=(m3/s)
Q_G=molar_volume(1)*Fvec(1);
%Liquid volume flow=sum((m3/mol)_i*(mol/s)_i)
Q_L=dot(molar_volume(2:6),Fvec(2:6));
end