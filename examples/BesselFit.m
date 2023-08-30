function a=  BesselFit(x,FitTarget, order)

[m,n]=size(x);
if n>1 
    x = x';
    FitTarget = FitTarget';
end

PlotFigure = 1;
% coefficient = LatticePotential*Erecoil/hbar/2;
ft = fittype(@(a,x)besselj(order,x*a).^2);
f = fit(x,FitTarget,ft,'Start',10^5);
if PlotFigure
plot(f,x,FitTarget)
xlabel('t/s')
ylabel('Normalized Proprtion')
end
a = f.a;
end