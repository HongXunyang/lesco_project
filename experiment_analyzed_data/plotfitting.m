function output = plotfitting(f)
    x = linspace(-100,150,300);
    A0 = f.A0;
    A1 = f.A1;
    A2 = f.A2;
    A3 = f.A3;
    A4 = f.A4;
    x0 = f.x0;
    x1 = f.x1;
    x2 = f.x2;
    x3 = f.x3;
    x4 = f.x4;
    res = f.res;
    res1 = f.res1;
    mu = f.mu;
    f0 = fVoigt(x,x0,A0,res,mu);
    f1 = fVoigt(x,x1, A1, res,mu);
    f2 = fVoigt(x,x2, A2, res,mu);
    f3 = fVoigt(x,x3, A3, res,mu);
    f4 = fGauss(x,x4,A4,res1);

    plot(x,f0,"LineWidth",1.5,'Color',[0.5, 0.5, 0.5, 0.5],'LineStyle','--');
    hold on
    area(x,f4, 'LineStyle','none', 'FaceAlpha',0.9, 'FaceColor',[200, 200, 204]/255);
    area(x,f3, 'LineStyle','none', 'FaceAlpha',0.9, 'FaceColor',[255, 153, 204]/255);
    area(x,f2, 'LineStyle','none', 'FaceAlpha',0.9, 'FaceColor',[204, 102, 255]/255);
    area(x,f1, 'LineStyle','none', 'FaceAlpha',0.9, 'FaceColor',[51, 102, 255]/255);
    plot(x,f0+f1+f2+f3+f4,"LineWidth",4,'Color',[255/255, 102/255, 255/255, 0.7]);
    output = 0;
end