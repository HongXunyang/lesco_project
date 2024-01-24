clear
clc
T21 = Fitting_fun('21', true);
T62 = Fitting_fun('62');
T104 = Fitting_fun('104');
T155 = Fitting_fun('155');

FIT.T21 = T21;
FIT.T62 = T62;
FIT.T104 = T104;
FIT.T155 = T155;

save('../results/FIT.mat','-struct','FIT');