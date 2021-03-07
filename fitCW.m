%% Fitting the CW penetration depth

%function fitCW()
close all
clear all
SAVE_AB = 1;
LOAD_AB = 0;

addpath('C:\OneDrivePolimi\OneDrive - Politecnico di Milano\Beta\Programs\MatlabTools')
%global ii mua mus AA BB;
%% try a fit of Z(rho) for different rho, mus=1 and fixed mua
% $$Z(\rho,\mu_a) = A(\mu_a)\rho^{B(\mu_a)}$$


%% GENERATE OR LOAD LIBRARY OF SIMULATIONS
if LOAD_AB == 0
    rho = [5:5:100]; %mm
    mua = [5e-3:1e-3:0.05]; %mm-1
    mus = [0.5:0.25:5.0]; %mm-1
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
    load 'Zrhomua'
end
nr=numel(rho);
ns=numel(mus);
na=numel(mua);

%% FIT Z with formulae
Z_calc=zeros(ns,na,nr); % calculated depth (mm)
Z_err=zeros(ns,na,nr); % calculated error on depth (mm)

% perform FIT
[A_fitresult, A_gof] = createFitNew(mua, mus, AA, 'poly33', 'AA')
[B_fitresult, B_gof] = createFitNew(mua, mus, BB, 'poly33', 'BB')

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
% ox=order(1); oy=order(2); or=order(3); oc=order(4);
% nx=P(ox).n; ny=P(oy).n; nr=P(or).n; nc=P(oc).n;
% x=P(ox).a; y=P(oy).a; r=P(or).a; c=P(oc).a;
% lx=P(ox).l; ly=P(oy).l; lr=P(or).l; lc=P(oc).l;
% ux=P(ox).u; uy=P(oy).u; ur=P(or).u; uc=P(oc).u;

%x=rho
%y=Z or Err
indMua=[1,3,6,10,20];

nc=numel(indMua);
nr=3;

figure('Name','ThisName','Position',[0 0 scrsz(3) scrsz(4)]);
levels=1:100;
subplot1(nr,nc, 'XTickL', 'Margin', 'YTickL', 'Margin','YScale','linear');

x_title={'Zmax (mm)',' err (mm)','err rel(%)'};

for ir=1:nr
    for ic=1:nc
        subplot1(ic+nc*(ir-1));
        x=rho;
        y=mus;
        
        if ir==1, zvalues=squeeze(ZRhoMuaMus(:,indMua(ic),:)); end
        if ir==2, zvalues=squeeze(Z_err(:,indMua(ic),:)); end
        if ir==3, zvalues=squeeze(100*Z_err(:,indMua(ic),:)./ZRhoMuaMus(:,indMua(ic),:)); end

        pcolor(x,y,zvalues); shading interp; colormap(pink); hold on;
        contour(x,y,zvalues,levels,'b','ShowText','on'); hold off;


        %contour(x,y,zvalues,[tresh_z 1E30],'k','LineWidth',2);
    %     if(min(x)<0.5*max(x)), minx=0; else minx=min(x); end
    %     if(min(y)<0.5*max(y)), miny=0; else miny=min(y); end
    %     xlim([minx max(x)]); ylim([miny max(y)]);
        if(ir==nr), xlabel('rho (mm)'); end
        %if(ic==1), ylabel(x_title(ir)); end
        if(ic==1), ylabel('Mus (mm-1)'); end
        if(ir==1), title(['Mua = ' num2str(mua(indMua(ic))) ' mm-1']); end
        hold off;
    end
end

%suptitle('Main Title');
% NameFig=['Results\' Type '_' P(ox).l '_' P(oy).l '_' P(or).l '_' P(oc).l '_' Tag];
% if SaveFig, save_figure(NameFig); end




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