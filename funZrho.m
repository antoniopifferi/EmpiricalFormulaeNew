function Zrho=funZrho(s0,rho,mua,musp,nin,nout,M)

s0=s0*1e-3;
rho=rho*1e-3;
mua=mua*1e3;
musp=musp*1e3;

A=Afunction(nout,nin);
D=1/(3*musp);
mueff=sqrt(mua/D);
ze=2*A*D;
zep=2*D;
zs=1/(musp);

% nest 2 lines are equivalent to Contini equation (45)
m=(-M:M);
Rs0rho=fRs0rho(s0,rho,mua,musp,nin,nout,M);
% Rs0rhom0=fRs0rho(s0,rho,mua,musp,nin,nout,0);

m(m==0)=[];
zp0=2*m*s0+2*m*ze+2*m*zep+zs;
zn0=2*m*s0+2*m*ze+2*m*zep-2*ze-zs;
zp00=2*m*ze+2*m*zep+zs;
zn00=2*m*ze+2*m*zep-2*ze-zs;

zpm0=zs;
znm0=-2*ze-zs;

for i=1:length(rho)
    F=exp(-mueff*sqrt(rho(i)^2+zp0.^2))./(2*m.*sqrt(rho(i)^2+zp0.^2));
    G=exp(-mueff*sqrt(rho(i)^2+zn0.^2))./(2*m.*sqrt(rho(i)^2+zn0.^2));
    H=exp(-mueff*sqrt(rho(i)^2+zp00.^2))./(2*m.*sqrt(rho(i)^2+zp00.^2));
    I=exp(-mueff*sqrt(rho(i)^2+zn00.^2))./(2*m.*sqrt(rho(i)^2+zn00.^2));
    
    F0=zpm0*mueff*exp(-mueff*sqrt(rho(i)^2+zpm0^2))/(rho(i)^2+zpm0^2);
    G0=zpm0*exp(-mueff*sqrt(rho(i)^2+zpm0^2))/(rho(i)^2+zpm0^2)^(3/2);
    H0=znm0*mueff*exp(-mueff*sqrt(rho(i)^2+znm0^2))/(rho(i)^2+znm0^2);
    I0=znm0*exp(-mueff*sqrt(rho(i)^2+znm0^2))/(rho(i)^2+znm0^2)^(3/2);
    
    Zrho(i)=s0+sum(F-G-H+I)/(4*pi*Rs0rho(i))-s0/(4*pi*Rs0rho(i))*(F0+G0-H0-I0);
end

Zrho=Zrho*1e3;

% figure(1);%clf;
% plot(rho*1e3,Zrho,'r')



end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rs0rho=fRs0rho(s0,rho,mua,musp,nin,nout,M)

A=Afunction(nout,nin);
D=1/(3*musp);
mueff=sqrt(mua/D);
ze=2*A*D;
zs=1/(musp);
m=(-M:M);

for i=1:length(rho)
    zp0=2*m*(s0+2*ze)+zs;
    zn0=2*m*(s0+2*ze)-2*ze-zs;
    F=(  mueff./(rho(i)^2+zp0.^2) + 1./(rho(i)^2+zp0.^2).^(3/2)  ).*zp0.*exp(-mueff*sqrt(rho(i)^2+zp0.^2));
    G=(  mueff./(rho(i)^2+zn0.^2) + 1./(rho(i)^2+zn0.^2).^(3/2)  ).*zn0.*exp(-mueff*sqrt(rho(i)^2+zn0.^2));
    Rs0rho(i)=(1/4/pi)*sum(F-G);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A=Afunction(n1,n2)
% Computes parameter A; equation (27)
%
% n1: external medium
% n2: diffusing medium

% page 4590
n=n2/n1;

if n>1,
    
    % equations (30)
    t1 = 4*(-1-n^2+6*n^3-10*n^4-3*n^5+2*n^6+6*n^7-3*n^8-(6*n^2+9*n^6)*(n^2-1)^(1/2))/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    t2 = (-8+28*n^2+35*n^3-140*n^4+98*n^5-13*n^7+13*n*(n^2-1)^3*(1-(1/n^2))^(1/2))/...
        (105*n^3*(n^2-1)^2);
    t3 = 2*n^3*(3+2*n^4)*log(...
        ((n-(1+n^2)^(1/2))*(2+n^2+2*(1+n^2)^(1/2))*(n^2+(n^4-1)^(1/2)))/...
        (n^2*(n+(1+n^2)^(1/2))*(-n^2+(n^4-1)^(1/2)))...
        )/...
        ((n^2-1)^2*(n^2+1)^(7/2));
    t4 = ( (1+6*n^4+n^8)*log((-1+n)/(1+n))+4*(n^2+n^6)*log((n^2*(1+n))/(n-1)) )/...
        ((n^2-1)^2*(n^2+1)^3);
    
    % equation (29)
    B = 1+(3/2)*( 2*(1-1/n^2)^(3/2)/3+t1+t2+( (1+6*n^4+n^8)*(1-(n^2-1)^(3/2)/n^3) )/( 3*(n^4-1)^2) +t3 );
    C = 1-( (2+2*n-3*n^2+7*n^3-15*n^4-19*n^5-7*n^6+3*n^7+3*n^8+3*n^9)/(3*n^2*(n-1)*(n+1)^2*(n^2+1)^2) )-t4;
    A = B/C;
    
elseif n==1,
    
    %page 4591
    A=1;
    
else
    
    % equations (28)
    r1 = (-4+n-4*n^2+25*n^3-40*n^4-6*n^5+8*n^6+30*n^7-12*n^8+n^9+n^11)/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    r2 = (2*n^3*(3+2*n^4))/((n^2-1)^2*(n^2+1)^(7/2))*...
        log( (n^2*(n-(1+n^2)^(1/2)))*(2+n^2+2*(1+n^2)^(1/2))/...
        (n+(1+n^2)^(1/2))/(-2+n^4-2*(1-n^4)^(1/2)) );
    r3 = (4*(1-n^2)^(1/2)*(1+12*n^4+n^8))/...
        (3*n*(n^2-1)^2*(n^2+1)^3);
    r4 = ( (1+6*n^4+n^8)*log((1-n)/(1+n))+4*(n^2+n^6)*log((1+n)/(n^2*(1-n)))  )/...
        ((n^2-1)^2*(n^2+1)^3);
    
    % equation (27)
    A = (1+(3/2)*(8*(1-n^2)^(3/2)/(105*n^3))-(((n-1)^2*(8+32*n+52*n^2+13*n^3))/(105*n^3*(1+n)^2)+r1+r2+r3) )/...
        (1-(-3+7*n+13*n^2+9*n^3-7*n^4+3*n^5+n^6+n^7)/(3*(n-1)*(n+1)^2*(n^2+1)^2)-r4);
    
end

end


