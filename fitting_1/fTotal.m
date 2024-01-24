function y = fTotal(x, x2, x3, x4, A1, A2, A3, A4, res, res1, mu, Delta, Gamma, T)
    yexcitations = fS(x,A1,Delta,Gamma, T);
    yPhonon2 = fPVoigt(x,x2, A2, res,mu);
    yPhonon3 = fPVoigt(x,x3, A3, res,mu);
    yBack = fGauss(x,x4,A4,res1);
    y = yexcitations+yPhonon2+yPhonon3 + yBack;
end