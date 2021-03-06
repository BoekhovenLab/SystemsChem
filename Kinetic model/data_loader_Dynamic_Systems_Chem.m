function [HPLCtime, HPLCEDC, HPLCdimer, HPLCtrimer, HPLCtetramer, HPLCpentamer, HPLChexamer, HPLCNacylureamonomer] = data_loader_Dynamic_Systems_Chem(Fuel,Precursor,steadystate, seed)
   
HPLCtime = 0;
HPLCEDC = 0;
HPLCNacylureamonomer = 0;
HPLCdimer = 0;
HPLCtrimer = 0;
HPLCtetramer = 0;
HPLCpentamer = 0;
HPLChexamer = 0;
    
if Fuel == 100 && Precursor == 300 
Data = xlsread('Dynamic_Systems_Chem.xlsx',1); % Excel file name, number of sheet
Data = sortrows(Data,1);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end 
    
if Fuel == 75 && Precursor == 300 
Data = xlsread('Dynamic_Systems_Chem.xlsx',2); % Excel file name, number of sheet
Data = sortrows(Data,2);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end
                  
if Fuel == 50 && Precursor == 300 
Data = xlsread('Dynamic_Systems_Chem.xlsx',3); % Excel file name, number of sheet
Data = sortrows(Data,4);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end
            
if Fuel == 30 && Precursor == 300 
Data = xlsread('Dynamic_Systems_Chem.xlsx',4); % Excel file name, number of sheet
Data = sortrows(Data,5);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end
             
if Fuel == 0 && Precursor == 300 && steadystate == 1
Data = xlsread('Dynamic_Systems_Chem.xlsx',5); % Excel file name, number of sheet
Data = sortrows(Data,6);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end

if Fuel == 100 &&  Precursor == 300 && steadystate == 1
Data = xlsread('Dynamic_Systems_Chem.xlsx',6); % Excel file name, number of sheet
Data = sortrows(Data,2);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';        
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end

if Fuel == 75 &&  Precursor == 300 && steadystate == 1
Data = xlsread('Dynamic_Systems_Chem.xlsx',7); % Excel file name, number of sheet
Data = sortrows(Data,7);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end
           
if Fuel == 30 &&  Precursor == 300 && seed == 1
Data = xlsread('Dynamic_Systems_Chem.xlsx',8); % Excel file name, number of sheet
Data = sortrows(Data,8);
HPLCtime = Data(:,1).';
HPLCEDC = Data(:,2).';
HPLCNacylureamonomer = Data(:,3).';
HPLCdimer = Data(:,4).';
HPLCtrimer =  Data(:,5).';
HPLCtetramer = Data(:,6).';
HPLCpentamer = Data(:,7).';
HPLChexamer = Data(:,8).';
end
         