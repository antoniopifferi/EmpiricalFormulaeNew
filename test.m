function test()

figure(1);clf;

t=0:1e-1:6;
musp=1;
nin=1.4;
nout=1.4;
rho=30;
s0=[20000 40 20 10]
for i = 1:numel(s0)
    [Z,t0]=Zmean(t,s0(i),musp,nin,nout,rho);
    plot(t,Z);
    hold on;
end
legend(num2str(s0'))