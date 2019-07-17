function f=bnm2limNumSchem%(A)%params/y must match contunlim & plotnullclines
%Case 1 only (same TSs in each area), limit x only - upper or upper and lower
%
a=8; b=50;
u=.8;
%Fixed points (roots of p9 in order, then on limit X=u):
stab=[1];
sad=[2,4];
uns=[3];
ustab=[1,3];
usad=[2];
%}

basins=0;
ub=0;%1-u;%u bottom - <.5, u>.5
v=1/a;
vb=1/a-v;
inc=.05;%Quiver increment
qscale=.6;%Quiver scale
%
%Symbolic substitution:
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);
%Evaluate:
ci=subs(c,[as,gs,b1s,b2s],[a,1,b,b]);
fp=roots(fliplr(ci)); fp=double(fp); fp=round(fp,4);
fp(imag(fp)~=0)=[];
fp(fp<0)=[]; fp(fp>1)=[];
%f=fp;%length(fp);
fpx=sortrows(fp);
fpz=arrayfun(@(X) xnull(X,a,1),fpx);

%
%For upper limit u:
fpuz=roots([-2*a^2*b,3*a*b,-a-b,1-u]); fpuz=double(fpuz); fpuz=round(fpuz,4);
fpuz(imag(fpz)~=0)=[];
zupper=(1-u)*(1-a*u+2*a*u^2);
fpuz(fpuz<0)=[]; fpuz(fpuz>zupper)=[];
%f=fp;%length(fp);
fpuz=sortrows(fpuz);
fpux=u*ones(length(fpuz),1);
%}
%
eps=.005;%0.001 for pretty plots
X=(ub:eps:u);
Y=[(vb:eps/a:v)];%/2),(1/a/2+.01:.01:1/a)]; %check with k in contumlim & at bottom
lx=length(X);
ly=length(Y);
A=zeros(lx,ly);
epsilon=.001;

%pm=(a-sqrt(a^2-4*a*v))/(2*a); pp=(a+sqrt(a^2-4*a*v))/(2*a);
%qm=(b-sqrt(b^2-4*b*k*u))/(2*b*k); qp=(b+sqrt(b^2-4*b*k*u))/(2*b*k);

xoverFcn=@(t,Y)evUN(t,Y,u,ub,v,vb); 
options=odeset('Events',xoverFcn);

if basins==1
tout=[];
yout=[];
teout=[];
yeout=[];
ieout=[];
tend=200;
for i=1:lx
    Xi=X(i);
    yveci=zeros(1,ly);
    parfor j=1:ly
        Yj=Y(j);
        tstart=0;
        y0=[Xi,Yj];
        while tstart<tend
        [t,y,te,ye,ie]=ode23(@(t,Y)Zcontunlim(t,Y,a,b,u,ub,v,vb),[tstart,tend],y0,options);
        %Accumulate output:  
        nt=length(t);
        tout=[tout;t(2:nt)];
        yout=[yout;y(2:nt,:)];
        %{
        %Mark events:
        teout=[teout;te];
        yeout=[yeout;ye];
        ieout=[ieout;ie];
        %}
        y0=y(end,:);
        %y0(y0<0)=0;
        %{
        if y0(1)>u
            y0(1)=u;
        elseif y0(1)<ub
            y0(1)=ub;
        end
        %}
        tstart=t(end);
        end
        %
        %Find closest f.p.:
        fpdiff=sqrt((y(end,1)-fpx).^2+(y(end,2)-fpz).^2);
        [~,nfp]=min(fpdiff);
        aij=a*fpz(nfp)-fpx(nfp);
        yveci(j)=aij;%A(i,j)=aij;
        %{
        xend=y(end,:);
        aij=a*xend(2)-xend(1);
        yveci(j)=aij;%A(i,j)=aij;
        %}
    end
    A(i,:)=yveci;
end
else
    xquiv=(ub:inc:u);
    yquiv=(vb:inc/a:v);
    U=zeros(length(xquiv),length(yquiv));
    V=U;
    for i=1:length(xquiv)
        xi=xquiv(i);
        for j=1:length(yquiv)
            yj=yquiv(j);
            arrow=Zcontunlim(0,[xi,yj],a,b,u,ub,v,vb);
            U(i,j)=arrow(1);
            V(i,j)=arrow(2);
        end
    end
    U=rot90(U);
    U=flipud(U);
    V=rot90(V);
    V=flipud(V);
    [xarrow,yarrow]=meshgrid(xquiv,yquiv);
end

xymin=.162; xymax=.838;
yymin=.265/a; yymax=.735/a;

epsilon=.002;
xlocs=[0,xymin,fpx(2:4)',u,xymax+5*epsilon,1];
xlocsPlot=[u,fpx(2),xymin,fpx(3),u,xymax,u,fpx(4)];
xlabs={'0','X_1^-','X_1^b','X_1^c','X_1^a','u','X_1^+','1'};
ylocs=[fpuz(1),fpz(2),yymin,fpz(3),fpuz(2)-epsilon,yymax-epsilon/2,fpuz(3)+epsilon/2,fpz(4)+epsilon,1/a];
ylocsPlot=[fpuz(1),fpz(2),yymin,fpz(3),fpuz(2),yymax,fpuz(3),fpz(4),1/a];
ylabs={'Y_1^{ul}','Y_1^b','Y_1^-','Y_1^c','Y_1^{um}','Y_1^+','Y_1^{uu}','Y_1^a','1/\alpha  '};
%
f=A;
fs=12; lw=2; ms=10;
figure
colormap(.5+.5*redblue)
hold on

for i=1:length(xlocsPlot)
    plot([xlocsPlot(i),xlocsPlot(i)],[0,ylocsPlot(i)],'k--','linewidth',1)
    plot([0,xlocsPlot(i)],[ylocsPlot(i),ylocsPlot(i)],'k--','linewidth',1)
end

xx=(0:.001:1);
yy=(0:.001/a:1/a);
[p1,p2]=Zplotnullclines(a,b,xx,yy);
plot(xx,p1,'color',[0 0 .5],'LineWidth',lw)
if u<1
    h1=plot(xx(xx>u),p1(xx>u),'color',[1,1,1],'LineWidth',lw);
    h1.Color(4)=.8;
end
if u<1
    plot([u,u],[0,1/a],'k:','linewidth',lw)
    [top,~]=Zplotnullclines(a,b,u,0);
    plot([u,u],[0,top],'b-','linewidth',lw,'color',[0 0 .5])
end
if ub>0
    plot([ub,ub],[0,1/a],'k--','linewidth',lw)
    [top,~]=Zplotnullclines(a,b,ub,0);
    plot([ub,ub],[top,1/a],'-','linewidth',lw,'color',[0 0 .5])
end
if v<1/a
    plot([0,1],[v,v],'k--','linewidth',lw)
end
if vb>0
    plot([0,1],[vb,vb],'k--','linewidth',lw)
end
plot(p2,yy,'color',[.5 0 0],'LineWidth',lw)
if u<1
    h2=plot(p2(p2>u),yy(p2>u),'color',[1,1,1],'LineWidth',lw);
    h2.Color(4)=.8;
end
%Repeat:
if u<1
    plot([u,u],[0,1/a],'k:','linewidth',lw)
    [top,~]=Zplotnullclines(a,b,u,0);
    plot([u,u],[0,top],'b-','linewidth',lw,'color',[0 0 .5])
end
%
%Fixed points (roots of p9 in order):
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
%Limit (extra ones):
%ustab=[];
if isempty(ustab)==0
    plot(fpux(ustab),fpuz(ustab),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw)%[.75,.75,.75]
end
%usad=[];
if isempty(usad)==0
    plot(fpux(usad),fpuz(usad),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw)
    plot(fpux(usad),fpuz(usad),'kx','markersize',ms,'markerfacecolor','k','linewidth',lw)
end
%
hold off
set(gca,'XTick',xlocs,'xticklabels',xlabs)
set(gca,'YTick',ylocs,'yticklabels',ylabs)
set(gca,'FontSize',fs)
xlabel('X_1','FontSize',20)
ylabel('Y_1','rot',0,'FontSize',20)
axis ([0,1,0,1/a])
box on
%f=x;
end

function [value,isterminal,direction]=evUN(t,Y,u,ub,v,vb)
l=1;%length(Y);
value=[Y(1)-u;ub-Y(1);Y(2)-v;vb-Y(2)];%Y-[u];
isterminal=ones(4,1);
direction=ones(4,1);
end

function f=Zcontunlim(t,Y,a,b,u,ub,v,vb)
x=Y(1); y=Y(2);
f1=x^2*(1-x)-x*y-(1-x)^2*x+(1-x)*(1/a-y);
f2=b*y^2*(1-a*y)-x*y-a*b*(1/a-y)^2*y+(1-x)*(1/a-y);
if x>=u && f1>0
    f1=0;
end
if x<=ub && f1<0
    f1=0;
end
if y>=v && f2>0
    f2=0;
end
if y<=vb && f2<0
    f2=0;
end
f=[f1;f2];
end

function [f,g]=Zplotnullclines(a,b,X,Y)
f=(1-X).*(1/a-X+2*X.^2);
g=(1-a*Y).*(1-b*Y+2*a*b*Y.^2);
end

function f=xnull(X,a,g)
f=-(1+g)*X^3+(1+2*g)*X^2-(g+1/a)*X+1/a;
end