function f=mapping2to1

al=12;

k=4*al/(al-3);%10;
be=k/32*al.^3;%100;
xx=(0:.01:1);
%eta=(0:.005:.5);
%eta2=eta.^2;
eta2=(0:.005:.25);

HA=@(aa,bb) k./(2*aa.*bb*(aa.^2-6*aa+12)-k);
HB=@(aa,bb) k./(k+(2*(2*aa.^2-aa.*bb+2*bb))./(2*aa.*bb*(aa.^2-6*aa+12)-1));

c3=32*al.^3*be;
c2=16*al^2*be*(3-al);
c1=2*al*be*(al^2-6*al+12);
c0=2*(2*al^2-al*be+2*be);

a=HA(al,be);
b=HB(al,be);

%a=5;
%b=40;

d3=1;
d2=-2;
d1=(1+a)/a;
d0=(1-b)/(a*b);

f1=c3*eta2.^3+c2*eta2.^2+c1*eta2+c0;
f2=d3*xx.^3+d2*xx.^2+d1*xx+d0;

fs=15; lw=2;
figure
plot(xx,f2,'-','linewidth',lw)
hold on
plot(4*eta2,f1,'--','linewidth',lw)
plot([0,1],[0,0],'k-','linewidth',1)
xlabel('X, 4\eta^2')
ylabel('p_6','rot',0)
set(gca,'fontsize',fs)
axis tight%([0,1,-.1,.1])
legend('1 area','2 area')
grid on
grid minor
box on