function [k0, k1, k2, k3, k1_2, k1_3, k1_4, k1_5, k2_2, k2_3, k2_4, k2_5, k4_2, k4_3, k4_4, k4_5, k4_6, k5_1, k5_2, k5_3, k5_4, k5_5, s4, s5, s6] = data_loader_k_values_Dynamic_Systems_Chem(model, Assembly, steadystate, flux); %this will load all the data.
%this will load all the data.
i = 0;
hr= 3600;

k0 = 0;
k1 = 0;
k2 = 0; 
k3 = 0;
   
k1_2 = 0;
k1_3 = 0;
k1_4 = 0;
k1_5 = 0;
    
k2_2 = 0;
k2_3 = 0;
k2_4 = 0;
k2_5 = 0;
  
k4_2 = 0;
k4_3 = 0;
k4_4 = 0;
k4_5 = 0;
k4_6 = 0;
    
k5_1 = 0;
k5_2 = 0;
k5_3 = 0;
k5_4 = 0;
k5_5 = 0;
   
    
s4 = 0;
s5 = 0;
s6 = 0;
  
 
if model == 'A'
%sets all rate constants
k0 = 1E-5;         % 1st order direct hydrolysis of EDC  
%Activation reactions (k1)
k1 =  0.03;         % 2nd order reaction of x1COOH + EDC => x1COOEDC 
k1_2 = 0.033;    % 2nd order reaction of x2COOOC + EDC => x2COOEDC 
k1_3 = 0.15;      % 2nd order reaction of x3COOOC + EDC => x3COOEDC 
k1_4 = 0.15;      % 2nd order reaction of x4COOOC + EDC => x4COOEDC 
k1_5 = 0.15;      % 2nd order reaction of x5COOOC + EDC => x5COOEDC 
%Anhdride fromation reaction (k2)
%from x1COOEDC
k2 = 1;               % 2nd order reaction of x1COOEDC + x1COOH/x2COOOC/x3COOOC/x4COOOC/x5COOOC 
k2_2 = 1.1;        % 2nd order reaction of x2COOEDC + x1COOH/x2COOOC/x3COOOC/x4COOOC
k2_3 = 5;           % 2nd order reaction of x3COOEDC + x1COOH/x2COOOC/x3COOOC
k2_4 = 5;           % 2nd order reaction of x4COOEDC + x1COOH/x2COOOC
k2_5 = 5;           % 2nd order reaction of x5COOEDC + x1COOH
%N-acylurea formation  (k3) 
k3 = 0.013;        %  1st order reaction of x1COOEDC => x1COO + EDU N-Acylurea
%Anhydride hydrolysis (k5)
k4_2 = 4E-3;      % hydrolysis of x2COOOC
k4_3 = 4E-3;      % hydrolysis of x3COOOC
if Assembly == 1
k4_4 = 0.04E-3; % hydrolysis of x4COOOC
k4_5 = 0.04E-3; % hydrolysis of x5COOOC
k4_6 = 0.04E-3; % hydrolysis of x6COOOC
else
k4_4 = 4E-3;     % hydrolysis of x4COOOC
k4_5 = 4E-3;     % hydrolysis of x5COOOC
k4_6 = 4E-3;     % hydrolysis of x6COOOC
end
%transacetlyation rate constant (k5)
k5_1 = 40E-3;    % transacylation of 1 with 2-6 
k5_2 = 44E-3;    % transacylation of 2 with 2-6 
k5_3 = 100E-3;  % transacylation of 3 with 2-6 
k5_4 = 100E-3;  % transacylation of 4 with 2-6 
k5_5 = 100E-3;  % transacylation of 5 with 2-5
if Assembly == 1;
s4 = 0.01;
s5 = 0.01;
s6 = 0.01;
else
s4 = 1;
s5 = 1;
s6 = 1;
end
end
