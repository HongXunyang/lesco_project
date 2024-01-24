function y = fGauss(x,x0,A,res)
    k = res/2;
    sigma = k/(sqrt(2)*sqrt(log(2)));
    y = A*exp(-(x-x0).^2/(2*sigma^2));
end

