function f=bnm2basinLim(a,b,u,stab,sad,uns,ustab,usad)%,A)%(A)%params/y must match contunlim & plotnullclines
%Case 1 only (same TSs in each area), limit x only - upper or upper and lower
%{
a=6; b=34;
u=1;
%Fixed points (roots of p9 in order, then on limit X=u):
stab=[1,5];
sad=[2,4];
uns=[3];
ustab=[];
usad=[];
%}
basins=1;
ub=0;%1-u;%u bottom - <.5, u>.5
v=1/a;
vb=1/a-v;
inc=.05;%Quiver increment
qscale=.6;%Quiver scale
eps=.05;%Basins - 0.001 for pretty plots
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
stabfpx=[fpx(stab);fpux(ustab)];
stabfpz=[fpz(stab);fpuz(ustab)];
%}
%
X=(ub:eps:u);
Y=[(vb:eps/a:v)];%/2),(1/a/2+.01:.01:1/a)]; %check with k in contumlim & at bottom
lx=length(X);
ly=length(Y);
epsilon=.001;

%pm=(a-sqrt(a^2-4*a*v))/(2*a); pp=(a+sqrt(a^2-4*a*v))/(2*a);
%qm=(b-sqrt(b^2-4*b*k*u))/(2*b*k); qp=(b+sqrt(b^2-4*b*k*u))/(2*b*k);

xoverFcn=@(t,Y)evUN(t,Y,u,ub,v,vb); 
options=odeset('Events',xoverFcn);

if basins==1
%
A=zeros(lx,ly);
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
        fpdiff=sqrt((y(end,1)-stabfpx).^2+(y(end,2)-stabfpz).^2);
        [~,nfp]=min(fpdiff);
        aij=a*stabfpz(nfp)-stabfpx(nfp);
        yveci(j)=aij;%A(i,j)=aij;
        %{
        xend=y(end,:);
        aij=a*xend(2)-xend(1);
        yveci(j)=aij;%A(i,j)=aij;
        %}
    end
    A(i,:)=yveci;
end
    f=A;
    %}
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
%A=round(A,3);
%{
uniquA=unique(reshape(A,size(A,1)*size(A,2),1));
lu=length(uniquA);
uniVec=(1:lu)';
for i=1:lu
    A(A==uniquA(i))=uniVec(i);
end
%}
%f=A;
fs=20; lw=2; ms=10;
figure
colormap(.5+.5*redblue)
if basins==1
    %contourf(X,Y,A',5,'LineStyle','none');%'LineWidth',2);%,'LineStyle','none');%,'LineWidth',2)%
    imagesc(X,Y,A')
    caxis([-1,1])
    set(gca,'YDir','normal')
elseif basins==0
    quiver(xarrow,yarrow,U,V,qscale,'k','linewidth',1);
end
hold on
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
if u<1
    if ub>0
        set(gca,'XTick',[0,ub,u,1],'xticklabels',{'0','u_-','u_+','1'})
    else
        set(gca,'XTick',[0,u,1],'xticklabels',{'0','u','1'})
    end
else
    set(gca,'XTick',[0,1])
end
if v<1/a
    if vb>0
        set(gca,'Ytick',[0,vb,v,1],'yticklabels',{'0','v_-','v_+','1/\alpha'})
    else
        set(gca,'Ytick',[0,v,1],'yticklabels',{'0','v','1/\alpha'})
    end
else
    set(gca,'Ytick',[0,1/a],'yticklabels',{'0','1/\alpha'})
end
set(gca,'FontSize',fs)
xlabel('X_1')
ylabel('Y_1','rot',0)
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

%{
if u-x(end,1)<epsilon
            if x(end,2)<qm%<epsilon
                A(i,j)=1;%Blue
            %else
                %A(i,j)=2;%Grey
            end
        elseif x(end,1)<pm%epsilon
            A(i,j)=-1;%Red
        end
        %{
        if x(end,1)<epsilon && v-x(end,2)<eps
            A(i,j)=-1;%Red
        elseif u-x(end,1)<epsilon && x(end,2)<epsilon
            A(i,j)=1;%Blue
        elseif u-x(end,1)<epsilon
            A(i,j)=2;%Grey
        end
        %}
%}
%plot([0,1],[0,0],'color',[.5 0 0],'LineWidth',lw)
%plot([0,0],[0,1/a],'color',[0 0 .5],'LineWidth',lw)
%
%plot(u,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(0,1/a,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%
%plot(xo,yo,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(xf,yf,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(xo2,yo2,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(pm,v,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(u,qm,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%
%plot(u,a*u*(1-u),'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(b*v*(1-k*v),v,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)