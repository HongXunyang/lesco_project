% Define your data
name = 'T21';
TEMP = load(['../results/',name,'.mat']);
FIT = load('../results/FIT.mat');
fs = FIT.(name).fs;
Data = TEMP.DATA;
max_e = 150;
min_e = -60;
data = Data(4);
x_base = data.en*1000;
x_base = x_base(x_base>min_e & x_base<max_e);
len_x = length(x_base(x_base>min_e & x_base<max_e));
len_q = length(Data);
Q = zeros(len_x,len_q);
E = zeros(len_x,len_q);
I = zeros(len_x,len_q);
for ii = 1:len_q
data = Data(ii);
en = data.en*1000;
cts = data.cts(en>min_e & en<max_e);
en = en(en>min_e & en< max_e);
new_y = interp1(en,cts,x_base) - fVoigt(x_base,fs{ii}.x0, fs{ii}.A0,fs{ii}.res,fs{ii}.mu);
E(:,ii) = x_base;
I(:,ii) = new_y;
Q(:,ii) = (I(:,ii))./I(:,ii)*data.qh;

end
T21.Q = Q; T21.E = E; T21.I = I;
EnergyMap_noElastic.T21 = T21;



% Define your data
name = 'T62';
TEMP = load([name,'.mat']);
FIT = load('FIT.mat');
fs = FIT.(name).fs;
Data = TEMP.DATA;
max_e = 150;
min_e = -60;
data = Data(4);
x_base = data.en*1000;
x_base = x_base(x_base>min_e & x_base<max_e);
len_x = length(x_base(x_base>min_e & x_base<max_e));
len_q = length(Data);
Q = zeros(len_x,len_q);
E = zeros(len_x,len_q);
I = zeros(len_x,len_q);
for ii = 1:len_q
data = Data(ii);
en = data.en*1000;
cts = data.cts(en>min_e & en<max_e);
en = en(en>min_e & en< max_e);
new_y = interp1(en,cts,x_base) - fVoigt(x_base,fs{ii}.x0, fs{ii}.A0,fs{ii}.res,fs{ii}.mu);
E(:,ii) = x_base;
I(:,ii) = new_y;
Q(:,ii) = (I(:,ii))./I(:,ii)*data.qh;

end
T62.Q = Q; T62.E = E; T62.I = I;
EnergyMap_noElastic.T62 = T62;



% Define your data
name = 'T104';
TEMP = load([name,'.mat']);
FIT = load('FIT.mat');
fs = FIT.(name).fs;
Data = TEMP.DATA;
max_e = 150;
min_e = -60;
data = Data(4);
x_base = data.en*1000;
x_base = x_base(x_base>min_e & x_base<max_e);
len_x = length(x_base(x_base>min_e & x_base<max_e));
len_q = length(Data);
Q = zeros(len_x,len_q);
E = zeros(len_x,len_q);
I = zeros(len_x,len_q);
for ii = 1:len_q
data = Data(ii);
en = data.en*1000;
cts = data.cts(en>min_e & en<max_e);
en = en(en>min_e & en< max_e);
new_y = interp1(en,cts,x_base) - fVoigt(x_base,fs{ii}.x0, fs{ii}.A0,fs{ii}.res,fs{ii}.mu);
E(:,ii) = x_base;
I(:,ii) = new_y;
Q(:,ii) = (I(:,ii))./I(:,ii)*data.qh;

end
T104.Q = Q; T104.E = E; T104.I = I;
EnergyMap_noElastic.T104 = T104;


% Define your data
name = 'T155';
TEMP = load([name,'.mat']);
FIT = load('FIT.mat');
fs = FIT.(name).fs;
Data = TEMP.DATA;
max_e = 150;
min_e = -60;
data = Data(4);
x_base = data.en*1000;
x_base = x_base(x_base>min_e & x_base<max_e);
len_x = length(x_base(x_base>min_e & x_base<max_e));
len_q = length(Data);
Q = zeros(len_x,len_q);
E = zeros(len_x,len_q);
I = zeros(len_x,len_q);
for ii = 1:len_q
data = Data(ii);
en = data.en*1000;
cts = data.cts(en>min_e & en<max_e);
en = en(en>min_e & en< max_e);
new_y = interp1(en,cts,x_base) - fVoigt(x_base,fs{ii}.x0, fs{ii}.A0,fs{ii}.res,fs{ii}.mu);
E(:,ii) = x_base;
I(:,ii) = new_y;
Q(:,ii) = (I(:,ii))./I(:,ii)*data.qh;

end
T155.Q = Q; T155.E = E; T155.I = I;
EnergyMap_noElastic.T155 = T155;


save('../results/EnergyMap_noElastic.mat','-struct','EnergyMap_noElastic');

