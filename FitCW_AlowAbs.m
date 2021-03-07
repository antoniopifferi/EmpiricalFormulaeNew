%% try to fit A coefficient
clear all
close all
load Zrhomua
x = mua;
mua_max = 0.002;

mask = mua<mua_max;
x = x(mask);
a = a(mask);

a0 = -1.672;%7e-4;%.002;
b0 = 10000;%5*1e-3 %1e-2;
c0 = 1.672;%-20;%1e-2;
d = 0.001;%0.005;
e = 0.85;
f = 0;%1e-2;
g = 1e-2;
FIT = 1;

%strfun='a0.*((x+d).^(-1.5)).*exp(-b0./(x+d) - c.*(x+d) -f.*(x+d))+e';
strfun = 'a0*exp(-b0*x)+c0';

y = eval(strfun);
plot(x,a,'.',x,y),%ylim([1e-3 1e3])

if FIT == 1
  
    % guess the mathematical expression
    strfun = 'a*exp(-b*x)+c';
    %strfun = 'a*x^b+c'
    
    [xData, yData] = prepareCurveData( x, a );
    
    % Set up fittype and options.
    ft = fittype( strfun, 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.StartPoint = [a0 b0 c0];
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts )
    
    % Plot fit with data.
    
    figure( 'Name', 'Fit A(mua)' );
    subplot(2,1,1)
    plot( fitresult, xData, yData );
    legend( 'A vs. mua', 'Fit A(\mu_a)', 'Location', 'NorthEast' );
    % Label axes
    xlabel \mu_a
    ylabel A
    grid on
    subplot(2,1,2)
    plot( fitresult, xData, yData, 'residuals' );
    drawnow;
end
