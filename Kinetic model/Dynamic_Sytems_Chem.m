clear;hr= 3600; HPLCtime = 0; HPLCEDC = 0; HPLCdimer = 0; HPLCtrimer = 0; HPLCtetramer = 0; HPLCpentamer = 0; HPLChexamer = 0; HPLCNacylureamonomer = 0; clf;  %clears some previously loaded data
 
%continuous fueling
z = 0;
x = 0;
timeup = 10* hr;
timestart = 0.33*hr;
flux = 100;                   %energy flux in mM/hr
interval = 0.01;            %period of the pulse in hr
pulse = flux * interval;  %amount of EDC added in mM 
 
steadystate = ;          %continuous fueling 0 = off; 1 = on.
seed = 0;                % seed on = 1 /off = 0;
Fuel = 50;               % fill in 30, 40, 50, 55, 75, 100
 

%Set time frames 
if steadystate == 1;           
t = 24.*hr;
elseif  steadystate == 0 && Fuel > 55;
t = 3*hr;
else
t = 1*hr;
end
%Set when assemblies are present 
if Fuel > 55; 
Assembly = 1;
else 
Assembly = 0;
end
 
 
[k0, k1, k2, k3, k1_2, k1_3, k1_4, k1_5, k2_2, k2_3, k2_4, k2_5, k4_2, k4_3, k4_4, k4_5, k4_6, k5_1, k5_2, k5_3, k5_4, k5_5, s4, s5, s6] = data_loader_k_values_Dynamic_Systems_Chem(Assembly, steadystate, flux, seed); 
%this will load all the data.
 
%model A: basic model 
%model B: updated model 
%This dataloader load all the optimized k-values depending on the chosen reaction conditions 
 
 
%these lines set initial Concentrations [=] M
EDC(1) = Fuel/1000;
 
%Assumption:
%INTERmolecular reaction both acid moieties can react with EDC. 
%INTRAmolecular reaction only one acid moiety reacts with EDC.
 
Precursor = 300;                  % set the concentration acid in mM
x1COOH(1) = (Precursor/1000);  
 
x1COOEDC(1) = 0;             %reactive monomer
x2COOOC(1) = 0;              %dimer
x2COOEDC(1) = 0;             %reactive dimer
 
x3COOOC(1) = 0;              %trimer
x3COOEDC(1) = 0;             %reactive trimer

if seed == 1;
x4COOOC(1) = 0.18/1000;      %tetramer
x4COOEDC(1) = 0;             %reactive tetramer

x5COOOC(1) = 0.14/1000;      %pentamer
x5COOEDC(1) = 0;             %reactive pentamer
 
x6COOOC(1) = 0.04/1000;      %hexamer
x6COOEDC(1) = 0;             %reactive hexamer
else
    
x4COOOC(1) = 0;              %tetramer
x4COOEDC(1) = 0;             %reactive tetramer

x5COOOC(1) = 0;              %pentamer
x5COOEDC(1) = 0;             %reactive pentamer
 
x6COOOC(1) = 0;              %hexamer
x6COOEDC(1) = 0;             %reactive hexamer
end
EDU(1)= 0;                   %waste
x1COO(1) = 0;                % N-Acylisourea monomer
 
%describes all the reactions 
for i=1:t;
   
    %nomenclature: rx_y_z
    %x refers to the type of reaction (activation, formation, deactivation, transacylation).
    
    %in activation reactions:
    %y refers to the species that reacts/is activated. 
    %z refers to the attacked species. 
    
    %in deactivation reactions:
    %y refers to the species that is deactivated. 
    %z refers to the smallest cleavage product).
    
    %in transacylation reactions:
    %y refers to the acid that reacts with an anhydride. 
    %z refers to the anhydride that is atacked by the acid. 
    %i refers to the different positions in the anhydride molecule that can be
    %attacked by the acid
 
    %Calculate all reaction reates
    r0(i) = k0*EDC(i);                                       %M/s EDC => EDU (direct hydrolysis)
   
    %Activation reactions (k1)
    r1_1(i) = k1*EDC(i)*x1COOH(i);            %x1COOH + EDC => x1COOEDC  activation of monomer
    r1_2(i) = k1_2*EDC(i)*x2COOOC(i);         %x2COOOC + EDC => x2COOEDC  activation of dimer
    r1_3(i) = k1_3*EDC(i)*x3COOOC(i);         %x3COOOC + EDC => x3COOEDC  activation of trimer
    r1_4(i) = k1_4*EDC(i)*x4COOOC(i);         %x4COOOC + EDC => x4COOEDC  activation of tetramer
    r1_5(i) = k1_5*EDC(i)*x5COOOC(i);         %x5COOOC + EDC => x5COOEDC  activation of pentamer
  
   %Anhydride formation reaction (k2)
 
    %from ACTIVATED MONOMER x1COOEDC
    r2_1_1(i) = k2*x1COOEDC(i)*x1COOH(i);       %x1COOH + x1COOEDC => x2COOOC + EDU dimer  formation
    r2_1_2(i) = k2_2*x1COOEDC(i)*x2COOOC(i);    %x2COOOC + x1COOEDC => x3COOOC + EDU trimer formation    via activated monomer    
    r2_1_3(i) = k2_3*x1COOEDC(i)*x3COOOC(i);    %x3COOOC + x1COOEDC => x4COOOC + EDU tetramer formation via activated monomer
    r2_1_4(i) = k2_4*x1COOEDC(i)*x4COOOC(i);    %x4COOOC + x1COOEDC => x5COOOC + EDU pentamer formation via activated monomer
    r2_1_5(i) = k2_5*x1COOEDC(i)*x5COOOC(i);    %x5COOOC + x1COOEDC => x6COOOC + EDU 
       
 
    %from ACTIVATED DIMER x2COOEDC
    r2_2_1(i) = k2*x2COOEDC(i)*x1COOH(i);        %x1COOH + x2COOEDC => x3COOOC + EDU
    r2_2_2(i) = k2_2*x2COOEDC(i)*x2COOOC(i);     %x2COOOC + x2COOEDC => x4COOOC + EDU
    r2_2_3(i) = k2_3*x2COOEDC(i)*x3COOOC(i);     %x3COOOC + x2COOEDC => x5COOOC + EDU  
    r2_2_4(i) = k2_4*x2COOEDC(i)*x4COOOC(i);     %x4COOOC + x2COOEDC => x6COOOC + EDU
    
    %from ACTIVATED TRIMER x3COOEDC
    r2_3_1(i) = k2*x3COOEDC(i)*x1COOH(i);       %x1COOH + x3COOEDC => x4COOOC + EDU 
    r2_3_2(i) = k2_2*x3COOEDC(i)*x2COOOC(i);    %x2COOOC + x3COOEDC => x5COOOC + EDU 
    r2_3_3(i) = k2_3*x3COOEDC(i)*x3COOOC(i);    %x3COOOC + x3COOEDC => x6COOOC + EDU 
    
    %from ACTIVATED TETRAMER x4COOEDC
    r2_4_1(i) = k2*x4COOEDC(i)*x1COOH(i);       %x1COOH+ x4COOEDC => x5COOOC + EDU 
    r2_4_2(i) = k2_2*x4COOEDC(i)*x2COOOC(i);    %x2COOOC+ x4COOEDC => x6COOOC + EDU
    
    %from ACTIVATED PENTAMER x5COOEDC
    r2_5_1(i) = k2*x5COOEDC(i)*x1COOH(i);       %x1COOH+ x5COOEDC => x6COOOC + EDU
     
    %N-acylurea formation (k3)
    r3(i) = k3*x1COOEDC(i);                     %x1COOEDC => x1COO + x1COOEDC 
 
    %Anhydride hydrolysis (k5)
    %DIMER x2COOOC
    r4_2_1(i) = k4_2*x2COOOC(i);                %x2COOOC => x1COO + x1COO 
      
    %TRIMER x3COOOC
    r4_3_1(i) = 2*k4_3*x3COOOC(i);              %x3COOOC => x2COOOC + x1COOH 
   
    %TETRAMER x4COOOC
    r4_4_1(i) = 2*k4_4*x4COOOC(i);              %x4COOOC => x3COOOC + x1COOH 
    r4_4_2(i) = k4_4*x4COOOC(i);                %x4COOOC => x2COOOC + x2COOOC 
   
    %PENTAMER x5COOOC
    r4_5_1(i) = 2*k4_5*x5COOOC(i);              %x5COOOC => x4COOOC + x1COOH  
    r4_5_2(i) = 2*k4_5*x5COOOC(i);              %x5COOOC => x3COOOC + x2COOOC 
   
    %HEXAMER x6COOOC
    r4_6_1(i) = 2*k4_6*x6COOOC(i);              %x6COOOC => x5COOOC + x1COOH 
    r4_6_2(i) = 2*k4_6*x6COOOC(i);              %x6COOOC => x4COOOC + x2COOOC 
    r4_6_3(i) = k4_6*x6COOOC(i);                %x6COOOC => x3COOOC + x3COOOC 
 
 
%Anhydride formation reaction (k5) - Transacylation  
 %from MONOMER x1COOH 
   r5_1_3(i) =  2*k5_1*x1COOH(i)*x3COOOC(i);          %x1COOH + x3COOOC => x2COOOC + x2COOOC 
   r5_1_4(i) =  s4*3*k5_1*x1COOH(i)*x4COOOC(i);       %x1COOH + x4COOOC => x3COOOC + x2COOOC 
      r5_1_5(i) =  s5*2*k5_1*x1COOH(i)*x5COOOC(i);    %x1COOH + x5COOOC => x4COOOC + x2COOOC 
      r5_1_6(i) =  s6*2*k5_1*x1COOH(i)*x6COOOC(i);    %x1COOH + x6COOOC => x5COOOC + x2COOOC 
      r5_1_5_ii(i) =  s5*2*k5_1*x1COOH(i)*x5COOOC(i); %x1COOH + x5COOOC => x3COOOC + x3COOOC
      r5_1_6_ii(i) =  s6*3*k5_1*x1COOH(i)*x6COOOC(i); %x1COOH + x6COOOC => x4COOOC + x3COOOC 
      
 %from DIMER x2COOOC
      r5_2_2(i) = k5_2*x2COOOC(i)*x2COOOC(i);         %x2COOOC + x2COOOC => x3COOOC + x1COOH 
      r5_2_3(i) = 2*k5_2*x2COOOC(i)*x3COOOC(i);       %x2COOOC + x3COOOC => x4COOOC + x1COOH 
      r5_2_4(i) = s4*2*k5_2*x2COOOC(i)*x4COOOC(i);    %x2COOOC + x4COOOC => x5COOOC + x1COOH
      r5_2_5(i) = s5*2*k5_2*x2COOOC(i)*x5COOOC(i);    %x2COOOC + x5COOOC => x6COOOC + x1COOH
      r5_2_4_ii(i) = s4*2*k5_2*x2COOOC(i)*x4COOOC(i); %x2COOOC + x4COOOC => x3COOOC + x3COOOC
      r5_2_5_ii(i) = s5*2*k5_2*x2COOOC(i)*x5COOOC(i); %x2COOOC + x5COOOC => x4COOOC + x3COOOC
      r5_2_6_ii(i) = s6*2*k5_2*x2COOOC(i)*x6COOOC(i); %x2COOOC + x6COOOC => x5COOOC + x3COOOC 
      r5_2_6_iii(i) = s6*2*k5_2*x2COOOC(i)*x6COOOC(i);%x2COOOC + x6COOOC => x4COOOC + x4COOOC 
      
%from TRIMER x3COOOC
      r5_3_2(i) = k5_3*x3COOOC(i)*x2COOOC(i);         %x3COOOC + x2COOOC => x4COOOC + x1COOH 
      r5_3_3(i) = 2*k5_3*x3COOOC(i)*x3COOOC(i);       %x3COOOC + x3COOOC => x5COOOC + x1COOH 
      r5_3_4(i) = s4*2*k5_3*x3COOOC(i)*x4COOOC(i);    %x3COOOC + x4COOOC => x6COOOC + x1COOH
      r5_3_3_i(i) = 2*k5_3*x3COOOC(i)*x3COOOC(i);     %x3COOOC + x3COOOC => x4COOOC + x2COOOC 
      r5_3_4_i(i) = s4*2*k5_3*x3COOOC(i)*x4COOOC(i);  %x3COOOC + x4COOOC => x5COOOC + x2COOOC
      r5_3_5_i(i) = s5*2*k5_3*x3COOOC(i)*x5COOOC(i);  %x3COOOC + x5COOOC => x6COOOC + x2COOOC
      
      r5_3_5_iii(i) = s5*2*k5_3*x3COOOC(i)*x5COOOC(i);%x3COOOC + x5COOOC => x4COOOC + x4COOOC
      r5_3_6_iii(i) = s6*2*k5_3*x3COOOC(i)*x6COOOC(i);%x3COOOC + x6COOOC => x5COOOC + x4COOOC
   
%from TETRAMER x4COOOC
      r5_4_2(i) = s4*k5_4*x4COOOC(i)*x2COOOC(i);      %x4COOOC + x2COOOC => x5COOOC + x1COOH 
      r5_4_3(i) = s4*2*k5_4*x4COOOC(i)*x3COOOC(i);    %x4COOOC + x3COOOC => x6COOOC + x1COOH 
      r5_4_3_i(i) = s4*2*k5_4*x4COOOC(i)*x3COOOC(i);  %x4COOOC + x3COOOC => x5COOOC + x2COOOC
      r5_4_4_i(i) = s4*2*k5_4*x4COOOC(i)*x4COOOC(i);  %x4COOOC + x4COOOC => x6COOOC + x2COOOC
      r5_4_4_ii(i) = s4*2*k5_4*x4COOOC(i)*x4COOOC(i); %x4COOOC + x4COOOC => x5COOOC + x3COOOC
      r5_4_5_ii(i) = s5*2*k5_4*x4COOOC(i)*x5COOOC(i); %x4COOOC + x5COOOC => x6COOOC + x3COOOC
      r5_4_6_ii(i) = s6*2*k5_4*x4COOOC(i)*x6COOOC(i); %x4COOOC + x6COOOC => x5COOOC + x5COOOC
       
%from PENTAMER x5COOOC 
      r5_5_2(i) = s5*k5_5*x5COOOC(i)*x2COOOC(i);      %x5COOOC + x2COOOC => x6COOOC + x1COOH 
      r5_5_3_i(i) = s5*2*k5_5*x5COOOC(i)*x3COOOC(i);  %x5COOOC + x3COOOC => x6COOOC + x2COOOC 
      r5_5_4_ii(i) = s5*2*k5_5*x5COOOC(i)*x4COOOC(i); %x5COOOC + x4COOOC => x6COOOC + x3COOOC
      r5_5_5_iii(i) = s5*2*k5_5*x5COOOC(i)*x5COOOC(i);%x5COOOC + x5COOOC => x6COOOC + x4COOOC
 
 %Monitoring the rate constants during the reaction cycle
   %MONOMER
   %Formation of x via hydrolysis
    r2_1_sum(i+1) = r4_2_1(i) + r4_3_1(i) + r4_4_1(i) + r4_5_1(i) + r4_6_1(i); 
   %Degradation of x via fuel
    r4_1_sum(i+1) = 1*(r2_1_1(i) + r2_2_1(i) + r2_3_1(i) + r2_4_1(i) + r2_5_1(i));
    %Transacylation
    %Formation of x via transacylation - form all x
    r5_1_a(i+1) = r5_2_2(i) + r5_2_3(i) + r5_2_4(i) + r5_2_5(i) + r5_3_2(i) + r5_3_3(i) + r5_3_4(i) + r5_4_2(i) + r5_4_3(i) + r5_5_2(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_1_d(i+1) =  1*(r5_1_3(i) + r5_1_4(i) + r5_1_5(i) + r5_1_6(i) + r5_1_5_ii(i) + r5_1_6_ii(i));
    
    %DIMER
    %Formation of x via fuel
    r2_2_sum(i+1) = r2_1_1(i);
    %Formation of x via hydrolysis
    r4_2_x_sum(i+1)= r4_3_1(i)+r4_4_2(i)+r4_5_2(i)+r4_6_2(i);
    %Degradation of x via hydrolysis 
    r4_2_sum(i+1) = 1*(r4_2_1(i));
    %Degradation of x via fuel
    r2_2_x_sum(i+1) = r1_2(i)+ r2_1_2(i)+r2_2_2(i)+r2_3_2(i)+r2_4_2(i);
    %Transacylation
    %Formation of x via transacylation - form all x!
    r5_2_a(i+1) = r5_1_3(i) + r5_1_4(i) + r5_1_5(i) + r5_1_6(i) + r5_3_3_i(i) + r5_3_4_i(i) + r5_3_5_i(i) + r5_4_3_i(i) + r5_4_4_i(i) + r5_5_3_i(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_2_d(i+1) = 1*(r5_2_2(i) + r5_2_3(i) + r5_2_4(i) + r5_2_5(i) + r5_2_4_ii(i) + r5_2_5_ii(i) + r5_2_6_ii(i) + r5_2_6_iii(i) + r5_3_2(i) + r5_4_2(i) + r5_5_2(i)) ;
    
    %TRIMER
    %Formation of x via fuel
    r2_3_sum(i+1) = r2_1_2(i) + r2_2_1(i);
    %Formation of x via hydrolysis
    r4_3_x_sum(i+1) =  r4_4_1(i)+r4_5_2(i)+r4_6_3(i);
    %Degradation of x via hydrolysis 
    r4_3_sum(i+1) =  1*(r4_3_1(i));
    %Degradation of x via fuel
    r2_3_x_sum(i+1) = r1_3(i)+r2_1_3(i)+r2_2_3(i)+ r2_3_3(i);
    %Transacylation
    %Formation of x via transacylation - form all x!
    r5_3_a(i+1) = r5_1_4(i) + r5_1_5_ii(i) + r5_1_6_ii(i) + r5_2_2(i) + r5_2_4_ii(i) + r5_2_5_ii(i) + r5_2_6_ii(i) + r5_4_4_ii(i) + r5_4_5_ii(i)+ r5_5_4_ii(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_3_d(i+1) = 1*(r5_1_3(i) + r5_2_3(i) + r5_3_2(i) + r5_3_3(i) + r5_3_4(i) + r5_3_3_i(i) + r5_3_4_i(i) + r5_3_5_i(i) + r5_3_5_iii(i) +  r5_3_6_iii(i) + r5_4_3(i) + r5_4_3_i(i)+ r5_5_3_i(i));
    
    %TETRAMER
    %Formation of x via fuel
    r2_4_sum(i+1) = r2_1_3(i) + r2_3_1(i) + r2_2_2(i);
    %Formation of x via hydrolysis
    r4_4_x_sum(i+1) =  r4_5_1(i)+r4_6_2(i);
    %Degradation of x via hydrolysis 
    r4_4_sum(i+1) = 1*(r4_4_1(i) + r4_4_2(i));
    %Degradation of x via fuel
    r2_4_x_sum(i+1) = r1_4(i)+r2_1_4(i)+r2_2_4(i);
    %Transacylation
    %Formation of x via transacylation - form all x!
    r5_4_a(i+1) = r5_1_5(i) + r5_1_6_ii(i) + r5_2_3(i) + r5_2_5_ii(i) + r5_2_6_iii(i) + r5_3_2(i) + r5_3_3_i(i) + r5_3_5_iii(i)  + r5_3_6_iii(i) + r5_5_5_iii(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_4_d(i+1) = 1*(r5_1_4(i) + r5_2_4(i) + r5_2_4_ii(i) + r5_3_4(i) + r5_3_4_i(i) + r5_4_2(i) + r5_4_3(i) + r5_4_3_i(i) + r5_4_4_i(i) + r5_4_4_ii(i) +  r5_4_5_ii(i) +  r5_4_6_ii(i) +  r5_5_4_ii(i));
   
    %PENTAMER
    %Formation of x via fuel
    r2_5_sum(i+1) = r2_1_4(i) + r2_4_1(i) + r2_2_3(i) + r2_3_2(i);
    %Formation of x via hydrolysis
    r4_5_x_sum(i+1) =  r4_6_1(i);
    %Degradation of x via hydrolysis 
    r4_5_sum(i+1) = 1*(r4_5_1(i) + r5_5_2(i));
    %Degradation of x via fuel
    r2_5_x_sum(i+1) = r1_5(i)+r2_1_5(i);
    %Transacylation
    %Formation of x via transacylation - form all x!
    r5_5_a(i+1) = r5_1_6(i) + r5_2_4(i) + r5_2_6_ii(i) + r5_3_3(i) + r5_3_4_i(i) + r5_3_6_iii(i) + r5_4_2(i) + r5_4_3_i(i) + r5_4_4_ii(i) + r5_4_6_ii(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_5_d(i+1) = 1*(r5_5_2(i) + r5_5_3_i(i) + r5_5_4_ii(i) + r5_5_5_iii(i) + r5_1_5(i) + r5_1_5_ii(i) + r5_2_5(i) + r5_2_5_ii(i) + r5_3_5_i(i) + r5_3_5_iii(i) + r5_4_5_ii(i));
   
    %HEXAMER
    %Formation of x via fuel
    r2_6_sum(i+1) = r2_1_5(i) + r2_5_1(i) + r2_2_4(i) + r2_4_2(i) + r2_3_3(i);
    %Degradation of x via hydrolysis 
    r4_6_sum(i+1) = 1*(r4_6_1(i) + r4_6_2(i) + r4_6_3(i));
    %Transacylation
    %Formation of x via transacylation - form all x!
    r5_6_a(i+1) = r5_2_5(i) + r5_3_4(i) + r5_3_5_i(i) + r5_4_3(i) + r5_4_4_i(i) + r5_4_5_ii(i) + r5_5_2(i) + r5_5_3_i(i) + r5_5_4_ii(i) + r5_5_5_iii(i);
    %Transacylation
    %Degradation of x via transacylation 
    r5_6_d(i+1) = 1*(r5_4_6_ii(i) + r5_3_6_iii(i) + r5_2_6_ii(i) + r5_2_6_iii(i) + r5_1_6(i)+ r5_1_6_ii(i)) ;
   
   %Calculation of the concentration at every second in the reaction cycle.
    %consumption of EDC when continuously fueling
    z = z+1;
    if z < interval*3600 && i < 0.33*hr && i < timeup && steadystate == 1 ;
     EDC(i+1) = EDC(i)-r0(i)-r1_1(i)-r1_2(i)-r1_3(i)-r1_4(i)-r1_5(i); %calculation of EDC concentration in M
    %calculation of EDC concentration in M
    z = 0;                                                   %reset timer
    elseif z > interval*3600 && 0.33*hr < i && i < timeup && steadystate == 1 ;
     EDC(i+1) = EDC(i)-r0(i)-r1_1(i)-r1_2(i)-r1_3(i)-r1_4(i)-r1_5(i) + (pulse/1000); 
    %calculation of EDC concentration in M
    z = 0;
    else %consumption of EDC when addition of a batch of EDC
    EDC(i+1) = EDC(i)-r0(i)-r1_1(i)-r1_2(i)-r1_3(i)-r1_4(i)-r1_5(i);
    end
    
    %consumption of monomer & formation of reactive monomer
    x1COOH(i+1) =  x1COOH(i)-r1_1(i)-r2_1_1(i)-r2_2_1(i)-r2_3_1(i)-r2_4_1(i)-r2_5_1(i)+2*r4_2_1(i)+r4_3_1(i)+r4_4_1(i)+r4_5_1(i)+r4_6_1(i)+r5_2_2(i)+r5_2_3(i)+r5_2_4(i)+r5_2_5(i)-r5_1_3(i)-r5_1_4(i)-r5_1_5(i)-r5_1_6(i)-r5_1_5_ii(i)-r5_1_6_ii(i)+r5_3_2(i)+r5_3_3(i)+r5_3_4(i)+r5_4_2(i)+r5_4_3(i)+r5_5_2(i);
    x1COOEDC(i+1) = x1COOEDC(i)+r1_1(i)-r2_1_1(i)-r2_1_2(i)-r2_1_3(i)-r2_1_4(i)-r2_1_5(i)-r3(i);   
     
    %formation of anhydride & reactive anhydrides
    %DIMER x2COOOC
    x2COOOC(i+1) =  x2COOOC(i)-r1_2(i)+r2_1_1(i)-r2_1_2(i)-r2_2_2(i)-r2_3_2(i)-r2_4_2(i)-r4_2_1(i)+r4_3_1(i)+2*r4_4_2(i)+r4_5_2(i)+r4_6_2(i)+2*r5_1_3(i)+r5_1_4(i)+r5_1_5(i)+r5_1_6(i)-2*r5_2_2(i)-r5_2_3(i)-r5_2_4(i)-r5_2_5(i)-r5_2_4_ii(i)-r5_2_5_ii(i)-r5_2_6_ii(i)-r5_2_6_iii(i)-r5_3_2(i)-r5_4_2(i)-r5_5_2(i)+r5_3_3_i(i)+r5_3_4_i(i)+r5_5_3_i(i)+r5_4_3_i(i)+r5_4_4_i(i)+r5_3_5_i(i);   
    x2COOEDC(i+1) = x2COOEDC(i)+r1_2(i)-r2_2_1(i)-r2_2_2(i)-r2_2_3(i)-r2_2_4(i); 
  
    %TRIMER x3COOOC
    x3COOOC(i+1) =  x3COOOC(i)-r1_3(i)+r2_2_1(i)+r2_1_2(i)-r2_1_3(i)-r2_2_3(i)-r2_3_3(i)-r4_3_1(i)+r4_4_1(i)+r4_5_2(i)+2*r4_6_3(i)+r5_1_4(i)-r5_1_3(i)+2*r5_1_5_ii(i)+r5_1_6_ii(i)+r5_2_2(i)-r5_2_3(i)+2*r5_2_4_ii(i)+r5_2_5_ii(i)+r5_2_6_ii(i)- r5_3_2(i)-2*r5_3_3(i)-r5_3_4(i)-2*r5_3_3_i(i)-r5_3_4_i(i)-r5_3_5_iii(i)-r5_3_6_iii(i)-r5_4_3(i)-r5_4_3_i(i)+r5_4_4_ii(i)+r5_4_5_ii(i)-r5_5_3_i(i)+r5_5_4_ii(i)-r5_3_5_i(i);
    x3COOEDC(i+1) = x3COOEDC(i)+r1_3(i)-r2_3_1(i)-r2_3_2(i)-r2_3_3(i); 
  
    %TETRAMER x4COOOC
    x4COOOC(i+1) = x4COOOC(i)-r1_4(i)+r2_1_3(i)+r2_2_2(i)+r2_3_1(i)-r2_1_4(i)-r2_2_4(i)-r4_4_1(i)-r4_4_2(i)+r4_5_1(i)+r4_6_2(i)+r5_1_5(i)-r5_1_4(i)+r5_2_3(i)-r5_2_4(i)-r5_2_4_ii(i)+r5_1_6_ii(i)+r5_2_5_ii(i)+2*r5_2_6_iii(i)+r5_3_2(i)-r5_3_4(i)+r5_3_3_i(i)-r5_3_4_i(i)+2*r5_3_5_iii(i)+r5_3_6_iii(i)-r5_4_2(i)-r5_4_3(i)-2*r5_4_4_i(i)-r5_4_3_i(i)-2*r5_4_4_ii(i)-r5_4_5_ii(i)-r5_4_6_ii(i)-r5_5_4_ii(i)+r5_5_5_iii(i);   
    x4COOEDC(i+1) = x4COOEDC(i)+r1_4(i)-r2_4_1(i)-r2_4_2(i);
  
    %PENTAMER x5COOOC
    x5COOOC(i+1) =  x5COOOC(i)-r1_5(i)+r2_1_4(i)+r2_2_3(i)+r2_3_2(i)+r2_4_1(i)-r2_1_5(i)-r4_5_1(i)-r4_5_2(i)+r4_6_1(i)+r5_1_6(i)-r5_1_5(i)-r5_1_5_ii(i)-r5_2_5(i)+r5_2_4(i)+r5_2_6_ii(i)-r5_2_5_ii(i)+r5_3_3(i)+r5_3_4_i(i)-r5_3_5_iii(i)+r5_3_6_iii(i)+2*r5_4_6_ii(i)+r5_4_2(i)+r5_4_3_i(i)+r5_4_4_ii(i)-r5_5_2(i)-r5_5_3_i(i)-r5_5_4_ii(i)-2*r5_5_5_iii(i)-r5_4_5_ii(i)-r5_3_5_i(i);  x5COOEDC(i+1) = x5COOEDC(i)+r1_5(i)-r2_5_1(i);
 
    %HEXAMER x6COOOC
    x6COOOC(i+1) =  x6COOOC(i)+r2_2_4(i)+r2_3_3(i)+r2_4_2(i)+r2_5_1(i)+r2_1_5(i)-r4_6_1(i)-r4_6_2(i)-r4_6_3(i)+r5_2_5(i)-r5_1_6(i)-r5_1_6_ii(i)-r5_2_6_ii(i)-r5_2_6_iii(i)+r5_3_4(i)-r5_3_6_iii(i)+r5_4_3(i)+r5_4_4_i(i)+r5_4_5_ii(i)-r5_4_6_ii(i)+r5_5_2(i)+r5_5_3_i(i)+r5_5_4_ii(i)+r5_5_5_iii(i)+r5_3_5_i(i);  
    x6COOEDC(i+1) = x6COOEDC(i);
  
    %formation of N-Acylurea 
    %N-MONOMER x1COO
    x1COO(i+1) =   x1COO(i)+r3(i);   
end
 
%load all data (if available)
[HPLCtime, HPLCEDC, HPLCdimer, HPLCtrimer, HPLCtetramer, HPLCpentamer, HPLChexamer, HPLCNacylureamonomer] = data_loader_Dynamic_Systems_Chem(Fuel,Precursor,steadystate, seed); %this will load all the data.
%Plots [EDC] vs time
    subplot(3,3,1)
    plot(HPLCtime,HPLCEDC,'x','color',[0.4940 0.1840 0.5560]);
    hold on
    plot((1:t)/60,EDC(1:t)*1000,'color',[0.4940 0.1840 0.5560]);  % plots EDC 
    hold on
    xlabel('Time (min)') % x-axis label
    ylabel('Concentration (mM)') % y-axis label
    title('EDC')
 
%Plots [monomer] vs time    
   subplot(3,3,2)
   hold on
   plot((1:t)/60,x1COOH(1:t)*1000,'-b');    % plots monomer 
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Monomer')
   hold on
  
 
%Plots [Anhydride] vs time    
   subplot(3,3,3)
   plot(HPLCtime,HPLCdimer,'x','color',[0.6350 0.0780 0.1840]);
   hold on
   plot((1:t)/60,x2COOOC(1:t)*1000,'color',[0.6350 0.0780 0.1840]);     % plots DIMER
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Dimer')
   hold on
   
   subplot(3,3,4)
   plot(HPLCtime,HPLCtrimer,'x','color',[0.9100    0.4100    0.1700]);
   hold on
   plot((1:t)/60,x3COOOC(1:t)*1000,'color',[0.9100    0.4100    0.1700]);  % plots TRIMER 
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Trimer')
   hold on
  
   subplot(3,3,5)
   plot(HPLCtime,HPLCtetramer,'x','color',[0 0.5 0]);
   hold on
   plot((1:t)/60,x4COOOC(1:t)*1000,'color',[0 0.5 0]);          % plots TETRAMER 
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Tetramer')
   hold on
  
   subplot(3,3,6)
   plot(HPLCtime,HPLCpentamer,'xm');
   hold on
   plot((1:t)/60,x5COOOC(1:t)*1000,'-m');    % plots pentamer 
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Pentamer')
   hold on
  
   subplot(3,3,7)
   plot(HPLCtime,HPLChexamer,'x','color',[0.25 0.25 0.25]);
   hold on
   plot((1:t)/60,x6COOOC(1:t)*1000,'color',[0.25 0.25 0.25]);           % plots HEXAMER
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('Hexamer')
   hold on
 
   %Plots [side product] vs time    
   subplot(3,3,8)
   plot(HPLCtime,HPLCNacylureamonomer,'xr');
   hold on
   plot((1:t)/60,x1COO(1:t)*1000,'-r');                 % plots N-acylisourea MONOMER 
   xlabel('Time (min)') % x-axis label
   ylabel('Concentration (mM)') % y-axis label
   title('N-acylurea monomer')
   hold on
   
   subplot(3,3,9)
   plot((1:t)/60,r2_3_sum(1:t)*1000,'r');    % plots formation of TRIMER via fuel
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
   
   subplot(3,3,9)
   plot((1:t)/60,r4_3_sum(1:t)*1000,'b');    % plots degradation of TRIMER via hydrolysis
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
   
   subplot(3,3,9)
   plot((1:t)/60,r5_3_a(1:t)*1000,'m');    % plots formation of TRIMER via transacylation
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
   
   subplot(3,3,9)
   plot((1:t)/60,r5_3_d(1:t)*1000,'c');    % plots degradation of TRIMER via transacylation
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
  
   subplot(3,3,9)
   plot((1:t)/60,r4_3_x_sum(1:t)*1000,'g');    % plots formation of TRIMER via hydrolysis
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
   
   subplot(3,3,9)
   plot((1:t)/60,r2_3_x_sum(1:t)*1000,'y');    % plots degradation of TRIMER via fuel
   xlabel('Time (min)') % x-axis label
   ylabel('Rate (mM^-s)') % y-axis label
   title('Activation Trimer')
   hold on
   
   
%sort out the data, selects every 60th datapoint (1min) and transposes it
%in order to work with it.
n = 10; 
r1_1 = round(transpose(1000000*r1_1(1 : n : end)),3); 
EDC = round(transpose(1000*EDC(1 : n : end)),2);
 
x2COOOC = round(transpose(1000*x2COOOC(1 : n : end)),5);
x3COOOC = round(transpose(1000*x3COOOC(1 : n : end)),5);
x4COOOC = round(transpose(1000*x4COOOC(1 : n : end)),5);
x5COOOC = round(transpose(1000*x5COOOC(1 : n : end)),5);
x6COOOC = round(transpose(1000*x6COOOC(1 : n : end)),5);
 
x1COO = round(transpose(1000*x1COO(1 : n : end)),5);
x1COOH = transpose(1000*x1COOH(1 : n : end));
EDU = transpose(1000*EDU(1 : n : end));
HPLCtime = transpose(HPLCtime);
HPLCNacylureamonomer = transpose (HPLCNacylureamonomer);
HPLCdimer = transpose(HPLCdimer);
HPLCtrimer = transpose(HPLCtrimer);
HPLCtetramer = transpose(HPLCtetramer);
HPLCpentamer = transpose(HPLCpentamer);
HPLChexamer = transpose(HPLChexamer);
 
r2_1_sum = round(transpose(1000*r2_1_sum(1 : n : end)),10);
r2_2_sum = round(transpose(1000*r2_2_sum(1 : n : end)),10);
r2_3_sum = round(transpose(1000*r2_3_sum(1 : n : end)),10);
r2_4_sum = round(transpose(1000*r2_4_sum(1 : n : end)),10);
r2_5_sum = round(transpose(1000*r2_5_sum(1 : n : end)),10);
r2_6_sum = round(transpose(1000*r2_6_sum(1 : n : end)),10);
 
r4_1_sum = round(transpose(1000*r4_1_sum(1 : n : end)),10);
r4_2_sum = round(transpose(1000*r4_2_sum(1 : n : end)),10);
r4_3_sum = round(transpose(1000*r4_3_sum(1 : n : end)),10);
r4_4_sum = round(transpose(1000*r4_4_sum(1 : n : end)),10);
r4_5_sum = round(transpose(1000*r4_5_sum(1 : n : end)),10);
r4_6_sum = round(transpose(1000*r4_6_sum(1 : n : end)),10);
 
r5_1_a = round(transpose(1000*r5_1_a(1 : n : end)),10);
r5_2_a = round(transpose(1000*r5_2_a(1 : n : end)),10);
r5_3_a = round(transpose(1000*r5_3_a(1 : n : end)),10);
r5_4_a = round(transpose(1000*r5_4_a(1 : n : end)),10);
r5_5_a = round(transpose(1000*r5_5_a(1 : n : end)),10);
r5_6_a = round(transpose(1000*r5_6_a(1 : n : end)),10);
 
r5_1_d = round(transpose(1000*r5_1_d(1 : n : end)),10);
r5_2_d = round(transpose(1000*r5_2_d(1 : n : end)),10);
r5_3_d = round(transpose(1000*r5_3_d(1 : n : end)),10);
r5_4_d = round(transpose(1000*r5_4_d(1 : n : end)),10);
r5_5_d = round(transpose(1000*r5_5_d(1 : n : end)),10);
r5_6_d = round(transpose(1000*r5_6_d(1 : n : end)),10);
 
 
r4_2_x_sum = round(transpose(1000*r4_2_x_sum(1 : n : end)),10);
r4_3_x_sum = round(transpose(1000*r4_3_x_sum(1 : n : end)),10);
r4_4_x_sum = round(transpose(1000*r4_4_x_sum(1 : n : end)),10);
r4_5_x_sum = round(transpose(1000*r4_5_x_sum(1 : n : end)),10);
 
 
r2_2_x_sum = round(transpose(1000*r2_2_x_sum(1 : n : end)),10);
r2_3_x_sum = round(transpose(1000*r2_3_x_sum(1 : n : end)),10);
r2_4_x_sum = round(transpose(1000*r2_4_x_sum(1 : n : end)),10);
r2_5_x_sum = round(transpose(1000*r2_5_x_sum(1 : n : end)),10);
 
 


