function y = fVoigt(x,x0,A,res,mu)
    y = mu*fLoren(x,x0,A,res) + (1-mu)*fGauss(x,x0,A,res);
end
