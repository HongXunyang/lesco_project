clear

temperatures = [21, 62, 104, 155];
Deltas = zeros(1,4); Deltas_err = zeros(1,4); 
Gammas = zeros(1,4); Gammas_err = zeros(1,4); 
As = zeros(1,4);
mus = [0.3,0.22,0.4,0.4];

figure;
for ii = 1:4
    temperature = temperatures(ii);
    T_sub = load(['T',num2str(temperature),'_sub.csv']);
    energy_min = -50;
    energy_max = 150;
    gammas = zeros(1,13);
    deltas = zeros(1,13);
    
    x = T_sub(:,8*2-1);
    y = T_sub(:,6*2)+T_sub(:,9*2)+T_sub(:,8*2)+T_sub(:,7*2);
    y = y(x>energy_min & x<energy_max);
    x = x(x>energy_min & x<energy_max);
    weights = ones(size(x));
    weights(x<0) = 0.3;
    ft = fittype('fTotal(x, x2, x3, x4, A1, A2, A3, A4, res, res1, mu, Delta, Gamma, T)');
    
    %ft(A1,A2,A3,A4,Delta,Gamma,T,mu,res,res1,x2,x3,x4)
    %                 A1  A2  A3  A4   Delta   Gamma   T             mu   res   res1  
    f = fit(x, y, ft,...
        'StartPoint', [1  1   1   1    30      20  temperature*0.086  mus(ii)  22.5   500 50 80 500 ],...
        'Lower',      [0  0   0   0    0       0   temperature*0.086  mus(ii)  20   100 40 50 100],...
        'Upper',      [100 10  10  10   100    100  temperature*0.086 mus(ii)  25  10000 80 110 1000], ...
        'Weights', weights);
    
    % Create subplot
    subplot(3, 2, ii);
    plotfitting(f);
    plot(x, y, '.', 'MarkerSize', 19.5, 'Color', [180, 70, 40]/255);
    title(['T = ', num2str(temperature), 'K']);

    % Set labels and ticks conditionally
    if ii > 2 % Bottom row
        xlabel('Energy (meV)');
    else
        set(gca, 'XTickLabel', []);
    end
    if mod(ii, 2) == 1 % Left column
        ylabel('Intensity (a.u.)');
    else
        set(gca, 'YTickLabel', []);
    end
    
    set(gca, 'FontSize', 14);
    conf = confint(f);
    Deltas(ii) = f.Delta; Deltas_err(ii) = f.Delta - conf(1,5);
    Gammas(ii) = f.Gamma; Gammas_err(ii) = f.Gamma - conf(1,6);
    As(ii) = f.A1;
end


subplot(3, 2, [5 6]); % Merge the bottom two cells
errorbar(temperatures, Deltas, Deltas_err, Deltas_err, '-o', 'LineWidth', 1.5, 'CapSize', 5); hold on
errorbar(temperatures, Gammas, Gammas_err, Gammas_err, '-o', 'LineWidth', 1.5, 'CapSize', 5);
legend(["Delta", "Gamma"]);
xlabel('Temperature (K)'); ylabel('Energy (mev)');
set(gca, 'FontSize', 14);

print(As)
    
    
    

