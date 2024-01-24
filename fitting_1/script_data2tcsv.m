clear
Results = load('../results/Result.mat');
energy_min = -100;
energy_max = 5;
Data = Results.Data;

num = 13;
len = 1300;

data_matrix = ones(len,2*num)*100000;

for temp = ["T21", "T62","T104", "T155"]
    data = Data.(temp);
    for ii = 1:num
        en = data(ii).en; en=en*1000;
        cts = data(ii).cts;
        data_matrix(1:length(en),2*ii-1) = en;
        data_matrix(1:length(cts),2*ii) = cts;
    end
    csvwrite([char(temp),'.csv'], data_matrix);
end

