clear 
Fit = load('../results/FIT.mat');
Simulation = load('simulation_result.mat');
Data.T21 = load('T21.mat').DATA;
Data.T62 = load('T62.mat').DATA;
Data.T104 = load('T104.mat').DATA;
Data.T155 = load('T155.mat').DATA;
Energy_no_elastic  = load('EnergyMap_noElastic.mat');
Energy_with_elastic  = load('EnergyMap_Elastic.mat');
Result.Data = Data;
Result.Simulation = Simulation;
Result.Fit_phonon = Fit;
Result.Energy_no_elastic = Energy_no_elastic;
Result.Energy_with_elastic = Energy_with_elastic;

save('Result.mat','-struct','Result');
