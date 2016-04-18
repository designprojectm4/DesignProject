function [visc_gas,visc_liq]= viscosity(T,P,Fvec)
%Finds Kinematic Viscosities for use in Re and Ga calculations
%Inputs:  
%       mu_H: Hydrogen Gas dynamic viscosity [Pa s]
%       mu_L=[mu_G;mu_A;mu_P;mu_W;mu_POH]: Liquid dynamic viscosity [Pa s]
%       mass density[rho_H;rho_G;rho_A;rho_P;rho_W;rho_POH]
%--------------------------------------------------------------------------
[mu_H,mu_L] = find_mu(T);%Calls dynamic viscosities
[~,mass_density] = density(T,P);
%--------------------------------------------------------------------------
%Unit Conversions
    function centistokes=SI2cSt(SI)%converts [m2 s-1] to Centistokes [cSt]
        centistokes=SI*1e6;
    end
    function SI=cSt2SI(centistokes)%converts Centistokes [cSt] to [m2 s-1]
        SI=centistokes*1e-6;
    end
%--------------------------------------------------------------------------
%HydrogenGas
mass_density_H=mass_density(1);%rho_H
visc_gas=mu_H/mass_density_H;%Hydrogen Kinematic Viscosity: [m2 s-1]
%--------------------------------------------------------------------------
%Liquid Mixture
%Calculations using VBN: Refutas (2000)
%https://neutrium.net/fluid_flow/estimating-the-viscosity-of-mixtures/

mass_density_L=mass_density(2:6);%returns[rho_G;rho_A;rho_P;rho_W;rho_POH]
kinematic_viscosities_L=bsxfun(@rdivide,mu_L,mass_density_L);
%returns a vector kinematic viscosities of pure liquids

vi=transpose(kinematic_viscosities_L);%Converts vector to array
vi=SI2cSt(vi);%convert to CentiStokes
VBN_i=14.534*log(log(vi+0.8))+10.975;
%Calculates VBN: Viscosity Blending Number
VBN_i=transpose(VBN_i);%Converts VBN array to vector

%Mass fraction calculations
Fvec_L=Fvec(2:6);%liquid molar flows [mol/s]
Fvec_L=Fvec_L/1000;%[kmol/s]
RFM=[92.09382;74.08;76.09;18.01528;60.09502];%[G;A;P;POH]
Fvec_mass=bsxfun(@times,Fvec_L,RFM);%mass flow vector
x_i=Fvec_mass/sum(Fvec_mass);%mass fraction vector

%Liquid mixture calculation
VBN_mixture=dot(x_i,VBN_i);%VBN of the mixture!
visc_liq=exp(exp((VBN_mixture-10.975)/14.534))-0.8;%in CentiStokes
visc_liq=cSt2SI(visc_liq);
%liquid mixture kinematic viscosity: [m2 s-1]
end

