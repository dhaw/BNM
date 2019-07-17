function f=HPCbasinsNonlin(whichplot)
basins=1;
%Parameters:
n=1;%Keep
k=9;%Alpha
%k=(length(x)-1)/(length(y)-1);
rx=1;
ry=20;
a=rx;
b=ry;
%
%Variables:
inc=.1;%.05;%Quiver increment
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
%HF=@(xx) a*(1-(xx/n).^pol);%RP
%HG=@(yy) b*(1-(k*yy/n).^pol);
%HF=@(xx) a/(n^deg)*(n-xx).^deg;%RQ
%HG=@(yy) b*k^deg/(n^deg)*(n/k-yy).^deg;
%HF=@(xx) polyapprox(n,xx,a,alpha,deg);
%HG=@(yy) polyapprox(n/k,yy,b,beta,deg);
%HF=@(xx) a/(1-exp(-expc))*(exp(-expc.*xx/n)-exp(-expc));
%HG=@(yy) b/(1-exp(-expc))*(exp(-expc.*yy*k/n)-exp(-expc));
%%
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
%%
%Calculate f.p.s:
v=(-1:.1:1)'; lv=length(v);
ics=[v,flipud(v)/k];
fps=zeros(lv,2);
xdot=@(xx,yy,k) xx.^2.*HF(xx)-xx.*yy-(1-xx).^2.*HF(1-xx)+(1-xx).*(1/k-yy);
ydot=@(xx,yy,k) yy.^2.*HG(yy)-xx.*yy-(1/k-yy).^2.*HG(1/k-yy)+(1-xx).*(1/k-yy);
deriv=@(t,xy,k) [xdot(xy(1),xy(2),k),ydot(xy(1),xy(2),k)]';
%{
[xmesh,ymesh]=meshgrid(x,y);
xdotgrid=xdot(xmesh,ymesh,k);
ydotgrid=ydot(xmesh,ymesh,k);
zgrid=abs(xdotgrid)-abs(ydotgrid);%1-(1-xdotgrid).*(1-ydotgrid);
c=contour(xmesh,ymesh,zgrid,[0,0]);
%}
%{
yofx=@(xx) xx.^2.*HF(xx)-(1-xx).^2.*HF(1-xx)-xx/k+1/k;
xofy=@(yy) k*yy.^2.*HG(yy)-k*(1/k-yy).^2.*HG(1/k-yy)-k*yy+1;
toSolve=@(xx) xofy(yofx(xx))-xx;%Intersection of nullclines
for i=1:length(ics);
    fun=@(xx) toSolve(xx);
    fpx(i,:)=fsolve(fun,ics(i,:));
end
%}
fpx=[0,.5,4];
fpx=round(fpx,4);
%fpx=unique(fpx);
fpy=HF(fpx);
A=zeros(lx,ly);
%%
if basins==1
    %A=zeros(lx,lz);
    epsilon=.001;
    for i=1:lx
        for j=1:ly
            [t,y]=ode23(@(t,Y)deriv(t,Y,k),[0,2300],[x(i),y(j)]);%,options);
            %Find closest f.p.:
            %fpdiff=sqrt((y(end,1)-fpx).^2+(y(end,2)-fpz).^2);
            %[~,nfp]=min(fpdiff);
            %A(i,j)=a*fpz(nfp)-fpx(nfp);
            yend=y(end,:);
            A(i,j)=a*yend(2)-yend(1);
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
f=A;
end