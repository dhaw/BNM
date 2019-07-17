function f=Zbif2D
%a=10; c1=(a+1)/a;%alpha
B=(.01:.01:100); lb=length(B);
alpha=1.9;
A=alpha*ones(1,lb);%1.9/6/7.1/10
%A=sqrt(B)+.05;
%A=(.01:.01:50);
%B=A.^2; lb=length(B);
X=zeros(9,lb);

%Convolution:
for i=1:lb
    b=B(i);
    a=A(i); c1=(a+1)/a;
    p1=[-2*a^2*b,3*a*b,-a-b,1];%g(y)
    p2=[-2,3,-c1,1/a];%f(x)
    %
    r=p1(1);
   for k=2:length(p1)
      r=conv(r,p2);
      r(end)=r(end)+p1(k);
   end
   %}
   %r=sym2poly(subs(poly2sym(p1),poly2sym(p2)));
   r(end-1)=r(end-1)-1;%-1: ** (above)
   X(:,i)=roots(r);
end
f=X;
%
X(imag(X)~=0)=NaN;%ox=ones(size(X));
b=repmat(B,9,1); a=repmat(A,9,1);
Y=-2*X.^3+3*X.^2-(a+1)./a.*X+1./a;
detJ=(-6*X.^2+6*X-(a+1)./a).*(-6*a.*b.*Y.^2+6*b.*Y-(a+b)./a)-1./a;
trJ=-6*X.^2+6*X-(a+1)./a-6*a.*b.*Y.^2+6*b.*Y-(a+b)./a;
Xsad=X; Xuns=X;
Xsad(detJ>0)=NaN;
X(detJ<0)=NaN; X(trJ>0)=NaN;
Xuns(detJ<0)=NaN; Xuns(trJ<0)=NaN;
%{
betam=(9*a-2*a^2-2*sqrt(a*(a-3)^3))/(4-a);
rm=roots([1,-2,c1,(1-betam)/a/betam]); Xm=max(real(rm));
betap=(9*a-2*a^2+2*sqrt(a*(a-3)^3))/(4-a);
rp=roots([1,-2,c1,(1-betap)/a/betap]); Xp=min(real(rp));
%}
X(imag(X)~=0)=NaN;
fs=20; ms=10; lw=2;
figure
hold on
%plot([betam,betam],[0,Xm],'k--','linewidth',lw)
%plot([betap,betap],[0,Xp],'k--','linewidth',lw)%Get rid if a>4

h1=plot(B,Xuns,'.','color',[1,0,0],'markersize',ms,'markerfacecolor',[1,0,0]);
h2=plot(B,X,'.','color',[0,1,0],'markersize',ms,'markerfacecolor',[0,1,0]);
h3=plot(B,Xsad,'k.','markersize',ms,'markerfacecolor','k');
%plot(B,X,'k.','markersize',ms,'markerfacecolor','k')
hold off
%set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam,betap],'xticklabels',{'0','\beta_-','\beta_+'})
%set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam],'xticklabels',{'0','\beta_-'})
set(gca,'fontsize',fs,'ytick',(0:.2:1),'yticklabel',{'0','','','','','1'})%,'xtick',[0:2:B(end)])
axis([0,B(end),0,1])
xlabel('\beta')
ylabel('X_1^e','rot',0)
l=legend([h1(1),h2(1),h3(1)],'Unstable node','Stable node','Saddle','location','NW');%[hleg,hobj,hout,mout]
grid on
%grid minor
box on
%hobj(1).Children.MarkerSize=20;
%}