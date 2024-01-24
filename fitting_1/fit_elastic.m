function f = fit_elastic(x,y)
    weights = ones(size(x));
    ft = fittype('fPVoigt(x,x0,A,res,mu)');
    % ft(A,mu,res,x0,x)
    f = fit(x, y, ft,...
        'StartPoint', [10  0.3  22.5  0  ],...
        'Lower',      [0   0    22.5  -5 ],...
        'Upper',      [100  1   22.5  5  ],...
        'Weights', weights);
end





