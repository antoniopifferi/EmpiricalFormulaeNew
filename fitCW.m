%% Fitting the CW penetration depth

%function fitCW()
close all
clear all
SAVE_AB = 1;
LOAD_AB = 0;
global ii mua mus AA BB;
%% try a fit of Z(rho) for different rho, mus=1 and fixed mua
% $$Z(\rho,\mu_a) = A(\mu_a)\rho^{B(\mu_a)}$$

if LOAD_AB == 0
    rho = [5:5:100]; %mm
    mua = [1e-3:1e-3:0.05]; %mm-1
    mus = [0.25:0.25:5.0]; %mm-1
%     rho = [5:5:100];
%     mua = [0:1e-3:0.05];
%     mus = [0.25:0.25:5.0];
    ZRhoMuaMus=zeros(numel(mus),numel(mua),numel(rho));
    for is = 1:numel(mus)
        mus_is = mus(is);
        nin = 1.4;
        nout = 1;
        M = 500;
        s0 = 1000;
        for ii = 1:numel(mua)
            Zrho(ii,:) = funZrho(s0,rho,mua(ii),mus_is,nin,nout,M);
            ZRhoMuaMus(is,ii,:)=Zrho(ii,:);
            result = createFit(rho,Zrho(ii,:),mod(ii,100)==0);
            a(ii) = result.a;
            b(ii) = result.b;
            AA(ii,is)=a(ii);
            BB(ii,is)=b(ii);
        end
    end
else
    load 'Zrhomua'
end

%% New Section
[A_fitresult, A_gof] = createFitNew(mua, mus, AA)
[B_fitresult, B_gof] = createFitNew(mua, mus, BB)

for is=1:numel(mus)
    for ia=1:numel(mua)
        for ir=1:numel(rho)
            A_calc=A_fitresult(mua(ia),mus(is));
            B_calc=B_fitresult(mua(ia),mus(is));
            Z_calc(is,ia,ir)=A_calc*rho(ir).^B_calc;
            err(is,ia,ir)=ZRhoMuaMus(is,ia,ir)-Z_calc(is,ia,ir);
        end
    end
end
ErrArray=err(:);
figure, plot(ErrArray);
for is=1:numel(mus)
    figure, surf(squeeze(err(is,4:end,:)))
end

for is=1:numel(mua)
    errTrunc=err(is:end,:,:);
    maxerrS(is)=max(abs(errTrunc(:)));
    miderrS(is)=median(abs(errTrunc(:)));
end

for ir=1:numel(rho)
    errTrunc=err(:,:,ir:end);
    maxerrR(ir)=max(abs(errTrunc(:)));
    miderrR(ir)=median(abs(errTrunc(:)));
end

for ia=1:numel(mua)
    errTrunc=err(:,ia:end,:);
    maxerrA(ia)=max(abs(errTrunc(:)));
    miderrA(ia)=median(abs(errTrunc(:)));
end

figure, semilogy(mua,maxerrA,mua,miderrA), xlabel('mua (mm-1)'), ylabel('max error (mm)'), grid on;
figure, semilogy(mus,maxerrS,mus,miderrS), xlabel('mus (mm-1)'), ylabel('max error (mm)'), grid on;
figure, semilogy(rho,maxerrR,rho,miderrR), xlabel('rho (mm)'), ylabel('max error (mm)'), grid on;
err_mean=mean(abs(ErrArray));
err_median=median(abs(ErrArray));
            
save All;


% %% New Fit
% na=numel(mua);
% nr=numel(rho);
% x=zeros(na*nr,2);
% for ia=1:na
%     for ir=1:nr
%         x(ir+ia*(nr-1),1)=rho(ir);
%         x(ir+ia*(nr-1),2)=mua(ir);
%     end
% end
% y=Zrho;
% size(x)
% p = polyfitn(x,y,3);
% 
% if exist('sympoly') == 2
%   polyn2sympoly(p)
% end
% if exist('sym') == 2
%   polyn2sym(p)
% end
% 
% % Evaluate on a grid and plot:
% [xg,yg]=meshgrid(0:.05:1);
% zg = polyvaln(p,[xg(:),yg(:)]);
% surf(xg,yg,reshape(zg,size(xg)))
% hold on
% plot3(x(:,1),x(:,2),y,'o')
% hold off
% 
% 
% 
%% Previous program
figure,plot(mua,[a',b']),legend({'A(\mu_a)','B(\mu_a)'},'FontSize',16)
xlabel \mu_a
if SAVE_AB == 1
    save('Zrhomua','Zrho','rho','mua','a','b')
end
drawnow;
%% try a fit for A(mua)
% $$A(\mu_a) = a\exp(-\mu_a b)+c$$
result_a = createFit_a(mua,a)
drawnow;
%% try a fit for B(mua)
% $$ B(\mu_a) = a_{2}\mu_a^{b_2}+c_2 $$
result_b = createFit_b(mua,b)
drawnow;
size(mua);
size(mus);
size(AA);
figure, surf(mua,mus,AA'), xlabel('mua (mm-1)'), ylabel('mus (mm-1)'), zlabel('AA');
figure, surf(mua,mus,BB'), xlabel('mua (mm-1)'), ylabel('mus (mm-1)'), zlabel('BB');


%% fitting
function [fitresult,gof] = createFit(rho,Zrho,verbosity)
% guess the mathematical expression
%strfun = 'a*(1-exp(-b*x))';
global ii mua
strfun = 'a*x.^b';

[xData, yData] = prepareCurveData( rho, Zrho);

% Set up fittype and options.
ft = fittype( strfun, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.TolFun = 1e-8;
opts.StartPoint = [0.0927335925457062 0.122034016208589];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
if verbosity == 1
figure( 'Name', 'Fit Z(\rho)' );
subplot(2,1,1)
plot( fitresult, xData, yData );
legend( 'Zrho vs. rho', 'fit Z(\rho|\mu_a)', 'Location', 'NorthEast' );
title(['\mu_a = ',num2str(mua(ii)),' mm^{-1}'])
% Label axes
xlabel rho
ylabel Zrho
grid on

subplot(2,1,2)
plot( fitresult, xData, yData, 'residuals' );

end
end
%% fitting
function [fitresult,gof] = createFit_a(mua, a)
% guess the mathematical expression
strfun = 'a*exp(-b*x)+c';
%strfun = 'a*x^b+c'

[xData, yData] = prepareCurveData( mua, a );

% Set up fittype and options.
ft = fittype( strfun, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.7 4.6 1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

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

function [fitresult,gof] = createFit_b(mua, b)
% guess the mathematical expression
%strfun = 'a*exp(-b*x)+c';
strfun = 'a*(x+d)^b+c';
%mua = mua(2:end); %to avoid zero
%b = b(2:end);
[xData, yData] = prepareCurveData( mua, b );

% Set up fittype and options.
ft = fittype( strfun, 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.TolFun = 1e-08;
opts.TolX = 1e-08;
opts.StartPoint = [0.012 -0.3642 0.5731 6.206e-5];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.

figure( 'Name', 'Fit B(mua)' );
subplot(2,1,1)
plot( fitresult, xData, yData );
legend( 'B vs. mua', 'Fit B(\mu_a)', 'Location', 'NorthEast' );
% Label axes
xlabel \mu_a
ylabel B
grid on
subplot(2,1,2)
plot( fitresult, xData, yData, 'residuals' );
drawnow;
end