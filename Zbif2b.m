function f=Zbif2b
%a=10; c1=(a+1)/a;%alpha
a=9;%[(.001:.001:15),(15.1:.1:30)];
g=1;%(.001:.001:5);
b1=40;%(1:.1:120);
b2=(1:.1:120);
param=b2;
lp=length(param); %****************
X=zeros(9,lp);
%
%Convolution:
for i=1:lp
    %
    %Change for correct parameter:
    b2=param(i); %g=b1./b2; %****************
    %??p1=[-a*(a*b1+b2),a*b1+2*b2,-(b2/a+a),1];%x=g(z)
    p1=[-a^2*(b1+b2),a*(b1+2*b2),-(a+b2),1];
    %??p2=[-(1+g),1+2*g,-(g+1/a),1/a];%z=f(x)
    p2=[-(1+g),(1+2*g),-(g+1/a),1/a];
    %
    r=p1(1);
   for k=2:length(p1)
      r=conv(r,p2);
      r(end)=r(end)+p1(k);
   end
   %
   %r=sym2poly(subs(poly2sym(p1),poly2sym(p2)));
   r(end-1)=r(end-1)-1;%-1: ** (above)
   X(:,i)=roots(r);
end
%}

%
X(imag(X)~=0)=NaN;
P=repmat(param,9,1);
%
%Change for correct parameter:
b2=P; %g=b1./b2;%param(j); b1=b2; %****************
Z=-(1+g).*X.^3+(1+2*g).*X.^2-(g+1./a).*X+1./a;
%detJ=zeros(9,lp); trJ=detJ;
x=X; z=Z;
%
detJ=(-3*(1+g).*x.^2+2*(1+2*g).*x-(g+1./a)).*(-3*a.*(b1+b2).*z.^2+2.*(b1+2*b2).*z-(b2./a+1))-1./a;
trJ=-3*(1+g).*x.^2+2*(1+2*g).*x-(g+1./a)+-3*a.*(b1+b2).*z.^2+2.*(b1+2*b2).*z-(b2./a+1);
%
Xsad=X; Xuns=X;
Xsad(detJ>0)=NaN;
X(detJ<0)=NaN; X(trJ>0)=NaN;
Xuns(detJ<0)=NaN; Xuns(trJ<=0)=NaN;

X(imag(X)~=0)=NaN;
fs=20; ms=10; lw=2;
figure
hold on
%plot([betam,betam],[0,Xm],'k--','linewidth',lw)
%plot([betap,betap],[0,Xp],'k--','linewidth',lw)%Get rid if a>4

%dotted=[3.2,12,20,27]; dotted=[dotted;dotted]; ld=length(dotted); ends=[zeros(1,ld);ones(1,ld)];
%plot(dotted,ends,'--','color',[.5,.5,.5],'linewidth',lw)

h1=plot(param,Xuns,'.','color',[1,0,0],'markersize',ms,'markerfacecolor',[1,0,0]);
h2=plot(param,X,'.','color',[0,1,0],'markersize',ms,'markerfacecolor',[0,1,0]);
h3=plot(param,Xsad,'k.','markersize',ms,'markerfacecolor','k');
paramval=56;
plot([paramval,paramval],[0,1],'k--','linewidth',lw)
%plot(B,X,'k.','markersize',ms,'markerfacecolor','k')
hold off
%set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam,betap],'xticklabels',{'0','\beta_-','\beta_+'})
%set(gca,'fontsize',fs,'ytick',[0,1],'xtick',[0,betam],'xticklabels',{'0','\beta_-'})
%set(gca,'fontsize',fs,'ytick',[0,1])%,'xtick',[0:2:B(end)])
set(gca,'fontsize',fs,'ytick',(0:.2:1),'yticklabel',{'0','','','','','1'})%,'xtick',[0:2:B(end)])
axis([0,param(end),0,1])
%Change for correct parameter:
xlabel('\beta_2')%****************
ylabel('X_1^e','rot',0)
legend([h1(1),h2(1),h3(1)],'Unstable node','Stable node','Saddle','location','NE')%'NW');northeastoutside
grid on
%grid minor
box on
%}
%%
%{
syms xs zs as gs b1s b2s
%p9(X):
%p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
p=subs(-as*(as*b1s+b2s)*zs^3+(as*b1s+2*b2s)*zs^2-(as+b2s/as)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);
%Jacobian:
%J=jacobian([xs^2*(1-xs)-xs*zs-gs*xs*(1-xs)^2+(1-xs)*(1/as-zs),b1s*zs^2*(1-as*zs)-xs*zs-b2s*zs*(1/a-zs)^2+(1-xs)*(1/as-zs)],[xs,zs]);
J=[2*(1+2*gs)*xs-3*(1+gs)*xs^2-gs-1/as,-1;-1/a,2*(b1s+2*b2s/as)*zs-3*(as*b1s+b2s)*zs^2-b2s/as^2-1];
%detJ=det(J);
dJ=(2*(1+2*gs)*xs-3*(1+gs)*xs^2-gs-1/as)*(2*(b1s+2*b2s/as)*zs-3*(as*b1s-b2s)*zs^2-b2s/as^2-1)-1/as;
tJ=2*(1+2*gs)*xs-3*(1+gs)*xs^2-gs-1/as+2*(b1s+2*b2s/as)*zs-3*(as*b1s-b2s)*zs^2-b2s/as^2-1;
%}
%%
%{
%Convolution:
for i=1:lp
    parami=param(i);
    %
    %Change for correct parameter:
    r=subs(c,[as,gs,b1s,b2s],[a,g,parami,parami]);
    X(:,i)=roots(fliplr(r));
end
f=X;
%}
%{
%for i=1:9
    %for j=1:lp
        x=X; z=Z;
        %
        %Change for correct parameter:
        b2=P; b1=b2;%param(j); b1=b2;
        %jac=subs(J,[xs,zs,as,gs,b1s,b2s],[X(i,j),Z(i,j),a,g,b1,b2]);
        %??detJ=(-3*(1+g)*x.^2+2*(1+2*g).*x-(g+1./a)).*(-3*(a*b1+b2).*z.^2+2*(b1+2*b2./a).*z-(b2./a^2+1))-1./a;
        detJ=(-3*(1+g).*x.^2+2*(1+2*g).*x-(g+1./a)).*(-3*(b1+b2).*z.^2+2.*(b1+2*b2).*z-(b2./a+1))-1./a;
        %detJ(i,j)=subs(dJ,[xs,zs,as,gs,b1s,b2s],[X(i,j),Z(i,j),a,g,b1,pj]);
        %detJ(i,j)=det(jac);
        %??trJ=-3*(1+g).*x.^2+2*(1+2*g).*x-(g+1./a)-3*(a.*b1+b2).*z.^2+2*(b1+2*b2./a).*z-(b2./a^2+1);
        trJ=-3*(1+g).*x.^2+2*(1+2*g).*x-(g+1./a)+-3*(b1+b2).*z.^2+2.*(b1+2*b2).*z-(b2./a+1);
        %trJ(i,j)=subs(tJ,[xs,zs,as,gs,b1s,b2s],[X(i,j),Z(i,j),a,g,b1,pj]);
        %trJ(i,j)=trace(jac);
    %end
%end
%}