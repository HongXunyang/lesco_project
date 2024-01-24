function y = fTotal(x, x0, x1, x2, x3, x4, A0,...
    A1, A2, A3, A4, res, res1, mu)

    yElastic = fVoigt(x,x0,A0,res,mu);
    yPhonon1 = fVoigt(x,x1, A1, res,mu);
    yPhonon2 = fVoigt(x,x2, A2, res,mu);
    yPhonon3 = fVoigt(x,x3, A3, res,mu);
    yBack = fGauss(x,x4,A4,res1);
    y = yElastic+yPhonon1+yPhonon2+yPhonon3 + yBack;
end