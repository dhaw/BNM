function f=plotTSetc
x=0:.01:1;
k0=4; k1=1;
y0=0:.001:1/k0;
y=0:.001:1/k1;
b=3.6;
%a1=2; a2=4; a3=8;
a1=1; a2=1; a3=1;

fs=20; lw=2; ms=10;
colx=[0,0,.5]; coly=[.5,0,0];
%
%Tolerance schedules
figure
hold on
plot([0,1/k0],[b,0],'-','linewidth',lw,'color',coly)
plot([0,1],[b,0],'-','linewidth',lw,'color',coly)
plot([1/k0,1],[0,0],'k--','linewidth',lw)
plot(1/k0,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(1/k1,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%txt1='P_1';
%text(1.5,7,txt1,'fontsize',fs)
set(gca,'FontSize',fs)
xlabel('Y')
ylabel('R_Y(Y)')%,'rot',0)
%grid on
%grid minor
box on
axis([0,1/k1,0,b])
set(gca,'xtick',[0,1/k0,1/k1],'xticklabels',{'0','1/k_0','1/k(t)'},'ytick',[0,b],'yticklabels',{'0','b'});
hold off
%}
%Nullclines
xofy0=b*y0.*(1-k0*y0);
xofy=b*y.*(1-k1*y);
%
xofy0=b*k0*y0.*(1-k0*y0);
xofy=b*k1*y.*(1-k1*y);
%}
yofx1=a1*x.*(1-x);
yofx2=a2*x.*(1-x);
yofx3=a3*x.*(1-x);

p1y=1/k0/2; p2y=1/k1/2;
p1x=b/4/k0; p2x=b/4/k1;
p1x=p2x;

figure
hold on
plot([0,1],1/k0*[1,1],'k-','linewidth',1)
h1=plot(x,yofx1,'k-','linewidth',lw,'color',colx);
h2=plot(x,yofx2,'k-','linewidth',lw,'color',colx);
h3=plot(x,yofx3,'k-','linewidth',lw,'color',colx);

h4=plot([0,1],[0,0],'k-','linewidth',lw,'color',colx);
plot([0,0],[0,1/k1],'k-','linewidth',lw,'color',coly);

h1.Color(4)=.2;
h2.Color(4)=.2;
h3.Color(4)=.2;
h4.Color(4)=.2;
plot(xofy0,y0,'k-','linewidth',lw,'color',coly)
plot(xofy,y,'k-','linewidth',lw,'color',coly)
plot(p1x,p1y,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)%,'color',coly);
plot(p2x,p2y,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)%,'color',coly);
plot([p1x,p2x],[p1y,p2y],'k--','linewidth',lw)%,'color',coly);
plot(0,1/k0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)%,'color',coly);
plot(0,1/k1,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)%,'color',coly);
plot([0,0],[1/k0,1/k1],'k--','linewidth',lw)%,'color',coly);
%txt1='P_1';
%text(1.5,7,txt1,'fontsize',fs)
set(gca,'FontSize',fs)
set(gca,'xtick',[0,1],'xticklabels',{'0','1'},'ytick',[0,1/k0,1/k1],'yticklabels',{'0','1/\alpha_0','1/\alpha(at)'});
xlabel('X')
ylabel('Z','rot',0)
grid on
grid minor
box on
axis([0,1,0,1/k1])
hold off