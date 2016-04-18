function H=reactionheat(T)
%--------------------------------------------------------------------------
%Conversion:
    function H_Joule=Cal2Joule(H_Cal)
       %Converts Calories to Joules
       H_Joule=H_Cal*4.184; 
    end

    function H_Joule=kCal2Joule(H_Cal)
        %Converts KiloCalories to Joules
        H_Joule=H_Cal*4184;
    end
%--------------------------------------------------------------------------
Tref=293;
%Formation Heat
T_form=[1;T*1e-2];
%Organics
Formations=[
    -124.03	-0.734%Glycerol
    -85.108	-0.931%Acetol
    -92.563	-0.903%Propylene Glycol
    -53.603	-0.938];%1-Propanol
H_form=Formations*T_form;%Cal/mol
H_form1=kCal2Joule(H_form);%J/mol
%Non Organics
T_form=[T;T^2;T^3;T^4;T^5];
Tref=[Tref;Tref^2;Tref^3;Tref^4;Tref^5];
IntC=[%integral of heat capacities wrt T
    1.6702,3.00E-04,1.00E-08,-1.00E-11,8.00E-16;%water
    13.759,0.0004,2.33E-07,-5.00E-11,4.00E-15];%hydrogen
H_form2_ref=[-241814;0];%standard temperature formation[water;hydrogen]
H_form2=H_form2_ref+IntC*T_form-IntC*Tref;
%List of formation enthalpies [J/mol]
HG=H_form1(1);
HA=H_form1(2);
HP=H_form1(3);
HPOH=H_form1(4);
HW=H_form2(1);
HH=H_form2(2);
%--------------------------------------------------------------------------
%Reaction Heat
H1=HA+HW-HG;%Heat of Reaction 1:G->A+W
H2=HP-HA-HH;%Heat of Reaction 2:A+H2->P
H3=HPOH+HW-HP-HH;%Heat of Reaction 3:P+H2->POH+H20
H=[H1;H2;H3];
end