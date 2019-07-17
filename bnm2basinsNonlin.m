function f=bnm2basinsNonlin(whichplot,A)
basins=1;
%Parameters:
n=1;%Keep
k=1;%Alpha
%k=(length(x)-1)/(length(y)-1);
rx=2;
ry=10;
a=rx;
b=ry;
fpx=[0,.5,1];
fpz=[1,.5,0]/k;
stab=[1,3];
sad=[2];
uns=[];
if whichplot==6
    fpx=[]; fpz=[];
    stab=[]; sad=[]; uns=[];
elseif whichplot==4
    %fpx=[]; fpz=[];
    stab=[2]; sad=[]; uns=[1,3];
elseif whichplot==2
    addfp=.122;
    fpx=[0,addfp,.5,1-addfp,1];
    fpz=[1,1-addfp,.5,addfp,0]/k;
    stab=[1,3,5];
    sad=[2,4];
    uns=[];
end
%
%Variables:
inc=.01;%.05;%Quiver increment
qscale=.6;%Quiver scale
x=(0:inc:1);
y=(0:inc/k:n/k);
lx=length(x);
ly=length(y);
x1=x;
y1=y;
%For functions:
pol=1;
pol2=2;
expc=4;
deg=3;%10000;
alpha=.75;
beta=.25;
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
%%
%
if whichplot==1
    HF=@(xx) a/(1-exp(-expc))*(exp(-expc.*xx/n)-exp(-expc));
    HG=@(yy) b*(1-(k*yy/n).^pol);
elseif whichplot==2
    HF=@(xx) a/(1-exp(-expc))*(exp(-expc.*xx/n)-exp(-expc));
    HG=@(yy) b/(1-exp(-expc))*(exp(-expc.*yy*k/n)-exp(-expc));
elseif whichplot==3
    HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
    HG=@(yy) b*(1-(k*yy/n).^pol);
elseif whichplot==4
    HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
    HG=@(yy) b*k^deg/(n^deg)*(n/k-yy).^deg;
elseif whichplot==5
    deg=2;
    HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
    HG=@(yy) b*(1-(k*yy/n).^pol);
elseif whichplot==6
    deg=2;
    HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
    HG=@(yy) b*k^deg/(n^deg)*(n/k-yy).^deg;
end
%}
%%
%Calculate f.p.s:
%{
v=(-1:.1:1)'; lv=length(v);
ics=[v,flipud(v)/k];
fps=zeros(lv,2);
%}
xdot=@(xx,yy,k) xx.^2.*HF(xx)-xx.*yy-(1-xx).^2.*HF(1-xx)+(1-xx).*(1/k-yy);
ydot=@(xx,yy,k) yy.^2.*HG(yy)-xx.*yy-(1/k-yy).^2.*HG(1/k-yy)+(1-xx).*(1/k-yy);
deriv=@(t,xy,k) [xdot(xy(1),xy(2),k),ydot(xy(1),xy(2),k)]';
%}
%Comment out if input A:
%{
A=zeros(lx,ly);
%%
if basins==1
    %A=zeros(lx,lz);
    epsilon=.001;
    for i=1:lx
        xi=x1(i);
        parfor j=1:ly
            [t,y]=ode23(@(t,Y)deriv(t,Y,k),[0,500],[xi,y1(j)]);%,options);
            yend=y(end,:);
            A(i,j)=k*yend(2)-yend(1);
        end
    end
else
    mx=max(x);
    my=max(y);
    u=zeros(length(x),length(y));
    v=u;
    for i=1:length(x)
        xi=x(i);
        for j=1:length(y)
            yj=y(j);
            u(i,j)=xi^2*HF(xi)-xi*yj-(1-xi)^2*HF(1-xi)+(1-xi)*(1/k-yj);
            v(i,j)=yj^2*HG(yj)-xi*yj-(1/k-yj)^2*HG(1/k-yj)+(1-xi)*(1/k-yj);
        end
    end
    u=rot90(u);
    u=flipud(u);
    v=rot90(v);
    v=flipud(v);
    [x,y]=meshgrid(x,y);
end
%}
f=A;
%%
fs=20; lw=2; ms=10;
figure
hold on
colormap(.5+.5*redblue)
if basins==1
    imagesc(x1,y1,A')
    caxis([-1,1])
    set(gca,'YDir','normal')
else
    quiver(x,y,u,v,qscale,'k','linewidth',1);
end
%axis([0,mx,0,my])
axis([0,n,0,n/k])
set(gca,'xtick',[0,n],'xticklabel',{'0','1'});
set(gca,'ytick',[0,n/k],'yticklabel',{'0','1/\alpha'});
xlabel('X_1')
ylabel('Y_1','rot',0)
set(gca,'FontSize',fs)

[xi,yi] = meshgrid(0:.001:n,0:.001:n/k);
HFX=arrayfun(HF,xi); HFX1=arrayfun(HF,1-xi);
HGY=arrayfun(HG,yi); HGY1=arrayfun(HG,1/k-yi);
g1=xi.^2.*HFX-xi.*yi-(1-xi).^2.*HFX1+(1-xi).*(1/k-yi);
g2=yi.^2.*HGY-xi.*yi-(1/k-yi).^2.*HGY1+(1-xi).*(1/k-yi);
g1(isreal(g1)==0)=NaN;
g2(isreal(g2)==0)=NaN;
contour(xi,yi,g1,[0 0],'color',[0 0 .5],'linewidth',lw);
if whichplot==6
    contour(xi,yi,g2,[0 0],'color',[.5 0 0],'linewidth',lw,'linestyle','--');
else
    contour(xi,yi,g2,[0 0],'color',[.5 0 0],'linewidth',lw);
end
%
%plot(fpx,fpy,'ko','markersize',ms,'markerfacecolor','k','linewidth',lw)
%Stable interior:
if isempty(stab)==0%stab=[] doesn't wotk if sad is a function arg
    plot(fpx(stab),fpz(stab),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw)
end
%Saddle interior:
if isempty(sad)==0
    plot(fpx(sad),fpz(sad),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw)
    plot(fpx(sad),fpz(sad),'kx','markersize',ms,'markerfacecolor','k','linewidth',lw)
end
%Unstable interior:
if isempty(uns)==0
    plot(fpx(uns),fpz(uns),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw)%[.75,.75,.75]
end
%plot(u,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(0,v,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
hold off
caxis([-1,1])
box on
end

function f=polyapprox(n,xx,a,alpha,deg)%a=b for Y
%ff=-a/alpha^deg*(xx/n)^deg+a;
%f=max([0,ff]);
if xx/n<alpha
    f=a;
else
    f=0;
end
end

%{
x=x1;
y=y1;
lx=length(x);
ly=length(y);
nx=zeros(lx,1);
ny=zeros(ly,1);
for i=1:lx
    nx(i)=x(i)*HF(x(i));
end
for i=1:ly
    ny(i)=y(i)*HG(y(i));
end
plot(x,nx,'color',[0 0 .5],'LineWidth',3)
plot(ny,y,'color',[.5 0 0],'LineWidth',3)
hold off
%}

%set(gca,'XTick',[0,1])
%set(gca,'YTick',[0,0.5,1])
%xlabel('X','FontSize',15)
%ylabel('Y','FontSize',15)

%{
[xi,yi] = meshgrid(-.1:.01:n,-.1:.01:n/k);
g1=xi.^2.*HF(xi)-xi.*yi;
g2=yi.^2.*HG(yi)-xi.*yi;
%g1(isreal(g1)==0)=NaN;
%g2(isreal(g2)==0)=NaN;
contour(xi,yi,g1,[0 0],'color',[0 0 .5],'linewidth',3);
contour(xi,yi,g2,[0 0],'color',[.5 0 0],'linewidth',3);
hold off
%}


%set(gca,'FontSize',22.5)%15)
%xlabel('X','FontSize',22.5)%15)
%ylabel('Y','FontSize',22.5)%15)