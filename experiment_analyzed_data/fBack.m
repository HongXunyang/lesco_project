function y = fBack(x,x4,A4,res1)
    y = -4*A4*x.*(x-x4)/(x4^2);
    y(y<0) = 0;
end