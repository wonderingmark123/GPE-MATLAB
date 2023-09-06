function a=  BesselFit(x,FitTarget, order , varargin)
if nargin == 4
    atry = varargin{1};
else 
    atry = 10^5;
end

[m,n]=size(x);
if n>1 
    x = x';
    FitTarget = FitTarget';
end

PlotFigure = 1;
% coefficient = LatticePotential*Erecoil/hbar/2;
ft = fittype(@(a,x)besselj(order,x*a).^2);
f = fit(x,FitTarget,ft,'Start',atry);
if PlotFigure
plot(f,x,FitTarget)
xlabel('t/s')
ylabel('Normalized Proprtion')
end
a = f.a;
end