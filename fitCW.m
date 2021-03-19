%% Fitting the CW penetration depth

close all
clear all
SAVE_AB = 1;
LOAD_AB = 0;
SAVE_FIG = 1;

%Parameters
indRho=[1,3,7,11,19];
MaxErr=5; % (mm) max displayed eaaror in plot
MaxPerc=10; % (%) max percentage relative error

addpath('C:\OneDrivePolimi\OneDrive - Politecnico di Milano\Beta\Programs\MatlabTools')
%global ii mua mus AA BB;
%% try a fit of Z(rho) for different rho, mus=1 and fixed mua
% $$Z(\rho,\mu_a) = A(\mu_a)\rho^{B(\mu_a)}$$


%% GENERATE OR LOAD LIBRARY OF SIMULATIONS
if LOAD_AB == 0
    rho = [5:2.5:50]; %mm
    mua = [3e-3:1e-3:0.05]; %mm-1
    mus = [0.5:0.25:5.0]; %mm-1
    ZRhoMuaMus=zeros(numel(mus),numel(mua),numel(rho));
    for is = 1:numel(mus)
        mus_is = mus(is);
        nin = 1.4;
        nout = 1;
        M = 500;
        s0 = 1000;
        for ia = 1:numel(mua)
            Zrho(ia,:) = funZrho(s0,rho,mua(ia),mus_is,nin,nout,M);
            ZRhoMuaMus(is,ia,:)=Zrho(ia,:);
            result = createFit(rho,Zrho(ia,:),mod(ia,100)==0);
            a(ia) = result.a;
            b(ia) = result.b;
            AA(ia,is)=a(ia);
            BB(ia,is)=b(ia);
        end
    end
else
    load 'All'
end
nr=numel(rho);
ns=numel(mus);
na=numel(mua);

RangeRho=[num2str(rho(1)) ',' num2str(rho(2)-rho(1)) ',' num2str(rho(nr))];
RangeMua=[num2str(mua(1)) ',' num2str(mua(2)-mua(1)) ',' num2str(mua(na))];
RangeMus=[num2str(mus(1)) ',' num2str(mus(2)-mus(1)) ',' num2str(mus(ns))];


%% ITERATE OVER POL ORDER
sRange=['Mua' RangeMua 'Mus' RangeMus 'Rho' RangeRho];
sMax=['Err' num2str(MaxErr) 'mmPerc' num2str(MaxPerc) '%'];
sDir=['Results_' sRange '_' sMax];
mkdir(sDir);

for iPa=1:3 % iterate over MUS POL Order
    for iPs=1:3 % iterate over MUS POL Order

%% FIT Z with formulae
Z_calc=zeros(ns,na,nr); % calculated depth (mm)
Z_err=zeros(ns,na,nr); % calculated error on depth (mm)

% perform FIT
polytype=['poly' num2str(iPa) num2str(iPs)];
[A_fitresult, A_gof] = createFitNew(mua, mus, AA, polytype, 'AA', sDir)
[B_fitresult, B_gof] = createFitNew(mua, mus, BB, polytype, 'BB', sDir)

% calc error
for is=1:ns
    for ia=1:na
        for ir=1:nr
            A_calc=A_fitresult(mua(ia),mus(is));
            B_calc=B_fitresult(mua(ia),mus(is));
            Z_calc(is,ia,ir)=A_calc*rho(ir).^B_calc;
            Z_err(is,ia,ir)=ZRhoMuaMus(is,ia,ir)-Z_calc(is,ia,ir);
        end
    end
end
save All;

figure, plot(abs(Z_err(:)));

%% Visualise ERROR
scrsz = get(0,'ScreenSize');

ncol=numel(indRho);
nrow=3;

FigError=figure('Name','ThisName','Position',[0 0 scrsz(3) scrsz(4)]);
levels=1:100;
subplot1(nrow,ncol, 'XTickL', 'Margin', 'YTickL', 'Margin','YScale','linear');

x_title={'Zmax (mm)',' err (mm)','err rel(%)'};

for ir=1:nrow
    for ic=1:ncol
        subplot1(ic+ncol*(ir-1));
        x=mua;
        y=mus;
        
        if ir==1, zvalues=squeeze(ZRhoMuaMus(:,:,indRho(ic))); end
        if ir==2, zvalues=squeeze(Z_err(:,:,indRho(ic))); end
        if ir==3, zvalues=squeeze(100*Z_err(:,:,indRho(ic))./ZRhoMuaMus(:,:,indRho(ic))); end

        pcolor(x,y,zvalues); shading interp; colormap(flipud(pink)); hold on;
        contour(x,y,zvalues,levels,'b','ShowText','on'); hold on;
        grid on;
        set(gca,'layer','top')

        if ir==2, caxis([0 MaxErr]); end
        if ir==3, caxis([0 MaxPerc]); end

        if(ir==nrow), xlabel('Mua (mm-1)'); end
        if((ic==1)&&(ir==1)), ylabel({'Zmax (mm)';'Mus (mm-1)'}); end
        if((ic==1)&&(ir==2)), ylabel({'err (mm)';'Mus (mm-1)'}); end
        if((ic==1)&&(ir==3)), ylabel({'rel err (%)';'Mus (mm-1)'}); end
        if(ir==1), title(['Rho = ' num2str(rho(indRho(ic))) ' mm']); end
        hold off;
    end
end

suptitle([sRange ' ' sMax ' ' polytype]);
NameFig=[sDir '\' 'Err_' polytype];
if SAVE_FIG, saveas(FigError,NameFig,'jpg'); end

    end % End iPs
end % End iPa




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