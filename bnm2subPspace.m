function f=bnm2subPspace%(A)
%Symbolic substitution:
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);

alpha=(1:.1:20); la=length(alpha);
beta=(2:.2:80); lb=length(beta);
A=zeros(la,lb);
%
for i=1:la
    a=alpha(i);
    for j=1:lb
        b=beta(j);
        cij=subs(c,[as,gs,b1s,b2s],[a,1,b,b]);
        fp=roots(fliplr(cij)); fp=double(fp); fp=round(fp,4);
        fp(imag(fp)~=0)=[];
        fp(fp<0)=[]; fp(fp>1)=[];
        A(i,j)=length(fp);
    end
end
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

figure
hold on
%cmap=[1,1,1;.5,.75,.75;.25,.5,.5];
cmap=[1,1,1;.4,.6,.6;0,.2,.2];
colormap(cmap);
contourf(alpha,beta,A','LineStyle','none');%'LineWidth',2);%,'LineStyle','none');%,'LineWidth',2)
%plot([2,2],[0,beta(end)],'k--','linewidth',lw)
plot(6,36,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%int1=0; plot([0,alpha(end)],[int1,2*alpha(end)+int1],'k--','linewidth',lw)
plot(alphaLine,betaLine,'k--','linewidth',lw)
plot(alphaLine2,betaLine2p,'k--','linewidth',lw)
plot(alphaLine2,betaLine2m,'k--','linewidth',lw)
%plot([0,alpha(end)],[16,16],'k--','linewidth',lw)
txt1='P';
text(2,9,txt1,'fontsize',20)
set(gca,'FontSize',fs)
xlabel('\alpha')
ylabel('\beta','rot',0)
grid on
grid minor
box on
ax = gca;
ax.Layer = 'top';
caxis([.5,9.5])
colorbar('ticks',[2,5,8],'ticklabels',{'3','5','9'})
axis([0,alpha(end),0,beta(end)])
hold off