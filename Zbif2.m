function f=Zbif2
a=5.9;%alpha
B=(1.01:.01:100); lb=length(B);
c1=(a+1)/a;
X=zeros(3,lb);
for i=1:lb
    b=B(i);
    X(:,i)=roots([1,-2,c1,(1-b)/a/b]);
end
X(imag(X)~=0)=NaN;
b=repmat(B,3,1);
Y=X.*(1-X);
detJ=(X.*(2-3*X)-Y).*(b.*Y.*(2-3*a*Y)-X)-X.*Y;
trJ=X.*(2-3*X)-Y+b.*Y.*(2-3*a*Y)-X;
Xsad=X;
Xsad(detJ>0)=NaN;
X(detJ<0)=NaN; X(trJ>0)=NaN;

betam=(9*a-2*a^2-2*sqrt(a*(a-3)^3))/(4-a);
rm=roots([1,-2,c1,(1-betam)/a/betam]); Xm=max(real(rm));
betap=(9*a-2*a^2+2*sqrt(a*(a-3)^3))/(4-a);
rp=roots([1,-2,c1,(1-betap)/a/betap]); Xp=min(real(rp));


fs=20; ms=10; lw=2;
figure
hold on
plot([betam,betam],[0,Xm],'k--','linewidth',lw)
%plot([betap,betap],[0,Xp],'k--','linewidth',lw)%Get rid if a>4
plot(B,Xsad,'k.','markersize',ms,'markerfacecolor','k')
plot(B,X,'.','color',[0,1,0],'markersize',ms,'markerfacecolor',[0,1,0])
hold off
%set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam,betap],'xticklabels',{'0','\beta_-','\beta_+'})
set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam],'xticklabels',{'0','\beta_-'})
axis([0,B(end),0,1])
xlabel('\beta')
ylabel('X','rot',0)