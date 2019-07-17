function f=pspaceAlphaBeta%(A)
alpha=(.05:.1:10.05); la=length(alpha);
beta=(2:.2:40); lb=length(beta);
A=zeros(la,lb);
[al,be]=meshgrid(alpha,beta);
al=al'; be=be';
A(be>2*al.^2./(al-2))=1;
A(al<2)=0;
m1=4./(8-al).*(9*al-al.^2+sqrt(al.*(al-6).^3));
m2=4./(8-al).*(9*al-al.^2-sqrt(al.*(al-6).^3));

B=zeros(la,lb);
be58=be; be58(al<5)=NaN; be58(al>8)=NaN;
%be5=be; be5(al<5)=NaN;
be8=be; be8(al<8)=NaN;
B(be58<m1)=.5;
B(be58>m2)=B(be58>m2)+.5;
B(be8>m2)=B(be8>m2)+1;
B(B>1)=1; B(B<1)=0;
A=A+B;
A(A==0)=NaN;
%

%}
f=A;
%Plot:
fs=20; lw=2; ms=10;
alphaLine=(2.01:.01:alpha(end));
betaLine=2*alphaLine.^2./(alphaLine-2);

alphaLine2=(6:.01:alpha(end));
qq=sqrt(alphaLine2.*(alphaLine2-6).^(3));
%alphaLine2a=(6:.01:8);%alpha(end));
%qqa=sqrt(alphaLine2a-6).^(3);
betaLine2p=(4./(alphaLine2-8)).*(alphaLine2.^2-9*alphaLine2+qq);
betaLine2m=(4./(alphaLine2-8)).*(alphaLine2.^2-9*alphaLine2-qq);

eps=.01; thr=5;%Max beta+
al=(3+eps:eps:alpha(end))';
lal=length(al);
[bet1,bet2]=arrayfun(@betfromal,al);

figure
hold on
%
cmap=[1,1,1;.67,.67,.67;.33,.33,.33];
colormap(cmap);
%contourf(alpha,beta,5*A','LineStyle','none');
%plot(3,9,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(6,36,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(alphaLine,betaLine,'k-','linewidth',lw)
%plot(alphaLine2,betaLine2p,'k-','linewidth',lw)
%plot(alphaLine2,betaLine2m,'k-','linewidth',lw)
%}
%
a=alpha; a(a<=2)=[]; a(a>=8)=[];
b=beta;%(beta>8);
%
%plot(alphaLine,betaLine,'k-','linewidth',lw)
%plot(alphaLine2,betaLine2p,'k-','linewidth',lw)
%plot(alphaLine2,betaLine2m,'k-','linewidth',lw)

%plot([0,alpha(end)],[1,1],'k-','linewidth',lw)
plot(al,bet1,'k--','linewidth',lw)
plot(al,bet2,'k--','linewidth',lw)

beta1=6;
plot([1,6],beta1*[1,1],'k-','linewidth',lw)
plot(6,beta1,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(1,beta1,'k<','markersize',ms,'markerfacecolor','k','linewidth',1)
beta2=16;
plot([4,8],beta2*[1,1],'k-','linewidth',lw)
plot(8,beta2,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(4,beta2,'k<','markersize',ms,'markerfacecolor','k','linewidth',1)
beta3=22;
plot([2,9],beta3*[1,1],'k-','linewidth',lw)
plot(9,beta3,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(2,beta3,'k<','markersize',ms,'markerfacecolor','k','linewidth',1)
beta4=34;
plot([1,7],beta4*[1,1],'k-','linewidth',lw)
plot(7,beta4,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(1,beta4,'k<','markersize',ms,'markerfacecolor','k','linewidth',1)

beta5=28;
plot([5,7],beta5*[1,1],'k-','linewidth',lw)
plot(7,beta5,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(5,beta5,'k<','markersize',ms,'markerfacecolor','k','linewidth',1)
%}
%txt1='P_1';
%text(1.5,7,txt1,'fontsize',fs)
%txt1='P_2';
%text(4.5,34,txt1,'fontsize',fs)
set(gca,'FontSize',fs)
xlabel('\alpha')
ylabel('\beta','rot',0)
grid on
%grid minor
box on
ax = gca;
ax.Layer = 'top';
%caxis([.5,9.5])
%colorbar('ticks',[2,5,8],'ticklabels',{'3','5','9'})
axis([0,alpha(end),0,beta(end)])
hold off
end

function f=toSolve(a,b)
f=b.^3-9*a.*b.^2+6*a.*(a+9).*b+16*a.^2;
end

function [f,g]=betfromal(al)
f=(9*al-2*al^2-2*sqrt(al*(al-3)^3))/(4-al);
g=(9*al-2*al^2+2*sqrt(al*(al-3)^3))/(4-al);
end