function f=Ztrajectories
%Parameters:
n=1;%Keep
k=3.1;
%k=(length(x)-1)/(length(y)-1);
rx=1;
ry=.2;
a=rx;
b=ry;
alpha=.75;
beta=.25;
%Variables:
x=(0:.1:1);
y=(0:.2:1/k);
x1=x;
y1=y;
%Functions:
pol=1;
pol2=2;
expc=4;
deg=10000;
%Functions:
%Use following for linear, with pol=1:
HF=@(xx) a*(1-(xx/n).^pol);%RP
HG=@(yy) b*(1-(k*yy/n).^pol);
%HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
%HG=@(yy) b*k^deg/(n^deg)*(n/k-yy).^deg;
%HF=@(xx) polyapprox(n,xx,a,alpha,deg);
%HG=@(yy) polyapprox(n/k,yy,b,beta,deg);
%HF=@(xx) a/(1-exp(-expc))*(exp(-expc.*xx/n)-exp(-expc));
%HG=@(yy) b/(1-exp(-expc))*(exp(-expc.*yy*k/n)-exp(-expc));

%HH=@(z) z*HF(z)*HG(z*HF(z));
HH=@(z) z-z*HF(z)*HG(z*HF(z));
ic=(0:.01:1); lic=length(ic);
z=zeros(lic,1);
for i=1:lic
    z(:,i)=fsolve(HH,ic(i));
end
%{
ali=a*k; betj=a*b;
xij=roots([1,-2,(1+ali)/ali,(1-betj)/(ali*betj)]); xij(imag(xij)~=0)=[];
xij=sortrows(xij); yij=xij.*(1-xij)*a;
if length(xij)==3
    xo1=xij(1); xo2=xij(3); xf=xij(2);
    yo1=yij(1); yo2=yij(3); yf=yij(2);
else
    xo1=xij(1); yo1=yij(1);
end
%}
%
z(z<0)=[]; z(z>1)=[];
z=round(z,4);
z=unique(z); z=sortrows(z,1);
zy=z.*arrayfun(HF,z);
%xo1=z(2); yo1=zy(2);%white
%xo2=z(4); yo2=zy(4);%white
%
%xf=z(3); yf=zy(3);%black
%}
mx=max(x);
my=max(y);
u=zeros(length(x),length(y));
v=u;
for i=1:length(x)
    xi=x(i);
    for j=1:length(y)
        yj=y(j);
        u(i,j)=xi^2*HF(xi)-xi*yj;
        v(i,j)=yj^2*HG(yj)-xi*yj;
    end
end
u=rot90(u);
u=flipud(u);
v=rot90(v);
v=flipud(v);
[x,y]=meshgrid(x,y);



fs=20; lw=2; ms=10;
figure
hold on
quiver(x,y,u,v,'k','linewidth',1.5);
axis([0,mx,0,my])
axis([0,n,0,n/k])
set(gca,'XTick',[0,n])
set(gca,'YTick',[0,n/k])
xlabel('X')
ylabel('Y','rot',0)
set(gca,'FontSize',fs)

[xi,yi] = meshgrid(0:.01:n,0:.01:n/k);
HFX=arrayfun(HF,xi);
HGY=arrayfun(HG,yi);
g1=xi.*HFX-yi;
g2=yi.*HGY-xi;
g1(isreal(g1)==0)=NaN;
g2(isreal(g2)==0)=NaN;
contour(xi,yi,g1,[0 0],'color',[0 0 .5],'linewidth',lw);
contour(xi,yi,g2,[0 0],'color',[.5 0 0],'linewidth',lw);

plot([0,n],[0,0],'color',[.5 0 0],'LineWidth',lw)
plot([0,0],[0,n/k],'color',[0 0 .5],'LineWidth',lw)

plot(1,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(0,1/k,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%
plot(xo1,yo1,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
plot(xf,yf,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(xo2,yo2,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
hold off
%}
end