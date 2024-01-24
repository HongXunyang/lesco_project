function FIT = Fitting_fun(name, ispic)

if nargin == 1
    ispic = false;
end



% T = 21     mu = 0.3
% T = 62     mu = 0.22
% T = 104    mu = 0.4
% T = 155    mu = 0.4-0.45

Muset.T21 = 0.3;
Muset.T62 = 0.22;
Muset.T104 = 0.4;
Muset.T155 = 0.4;
mu = Muset.(['T',name]);


Data = load(['T',name, '.mat']);
Data = Data.DATA;

amp0 = []; lo0 = []; up0 = [];
amp1 = []; lo1 = []; up1 = []; 
amp2 = []; lo2 = []; up2 = []; 
amp3 = []; lo3 = []; up3 = [];
en0 = []; enl0 = []; enu0 = [];
en1 = []; enl1 = []; enu1 = [];
en2 = []; enl2 = []; enu2 = [];
en3 = []; enl3 = []; enu3 = [];
Hset= [];
ndim = length(Data);
fit_set = cell(1,ndim);

% remove the first momentum and the last momentum
for ii  = 1:ndim
data = Data(ii);
energy = data.en;
intensity = data.cts;
x = 1000*energy((energy>-0.1) & (energy<0.15));
y = intensity((energy>-0.1) & (energy<0.15));
weights = ones(size(x));
weights(x>-5 & x<15) = 2;
Hset(ii) = data.qh;
ft = fittype('fTotal(x, x0, x1, x2, x3, x4, A0, A1, A2, A3, A4, res, res1, mu)');

%                 A0,  A1,  A2,  A3,  A4,  mu,  res,  res1,  x0,  x1,  x2,  x3,  x4
f = fit(x, y, ft,...
    'StartPoint', [2  0.3   0.3  0.3 0.1,  mu    22.5   500     0    30    50   80   500],...
    'Lower',      [0  0     0    0   0     mu    22.5   200     -30  0    0     0   300],...
    'Upper',      [80 10    10   10  1     mu    22.5  10000   30  50    80   110  10000],...
    'Weights', weights);
fit_set{ii} = f;
conf = confint(f); % error bar
amp0(ii) = f.A0; lo0(ii) = f.A0 - conf(1,1); up0(ii) = conf(2,1)-f.A0;
amp1(ii) = f.A1; lo1(ii) = f.A1 - conf(1,2); up1(ii) = conf(2,2)-f.A1;
amp2(ii) = f.A2; lo2(ii) = f.A2 - conf(1,3); up2(ii) = conf(2,3)-f.A2;
amp3(ii) = f.A3; lo3(ii) = f.A3-conf(1,4); up3(ii) = conf(2,4)-f.A3;

en0(ii) = f.x0; enl0(ii) = f.x0-conf(1,9); enu0(ii) = conf(2,9)-f.x0;
en1(ii) = f.x1; enl1(ii) = f.x1-conf(1,10); enu1(ii) = conf(2,10)-f.x1;
en2(ii) = f.x2; enl2(ii) = f.x2-conf(1,11); enu2(ii) = conf(2,11)-f.x2;
en3(ii) = f.x3; enl3(ii) = f.x3-conf(1,12); enu3(ii) = conf(2,12)-f.x3;

if ispic
    figure()
    plotfitting(f);
    plot(x,y,'.','MarkerSize', 13);
    title(['h = ', num2str(data.qh)])
    xlabel("Energy Loss (meV)");
    xlim([-100,150]);
    set(gca, 'FontSize', 14)
    legend("Elastic","Background","Breathing phonon","Phonon 50","Phonon 30","Fitting curve", "data")
end
end


Phonon1.amp = amp1; Phonon1.lo = lo1; Phonon1.up = up1;
Phonon2.amp = amp2; Phonon2.lo = lo2; Phonon2.up = up2;
Phonon3.amp = amp3; Phonon3.lo = lo3; Phonon3.up = up3;
Elastic.amp = amp0; Elastic.lo = lo0; Elastic.up = up0;

Phonon1.en = en1; Phonon1.enl = enl1; Phonon1.enu = enu1;
Phonon2.en = en2; Phonon2.enl = enl2; Phonon2.enu = enu2;
Phonon3.en = en3; Phonon3.enl = enl3; Phonon3.enu = enu3;
Elastic.en = en0; Elastic.enl = enl0; Elastic.enu = enu0;


FIT.phonon1 = Phonon1;
FIT.phonon2 = Phonon2;
FIT.phonon3 = Phonon3;
FIT.elastic = Elastic;
FIT.h = Hset;
FIT.fs = fit_set;


end