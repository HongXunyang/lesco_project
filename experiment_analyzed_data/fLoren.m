function y = fLoren(x,x0,A,res)
    y = A*(res/2)^2./((x-x0).^2 + (res/2)^2);
end