function f=bnm2basinDiffTS%(A)%a=alpha, b=beta
%
basins=1;%>0, not 1 if input basins
a=9; g=1; b1=40; b2=56;
%Symbolic substitution:
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);
%Evaluate:
ci=subs(c,[as,gs,b1s,b2s],[a,g,b1,b2]);
fp=roots(fliplr(ci)); fp=double(fp); fp=round(fp,4);
fp(imag(fp)~=0)=[];
fp(fp<0)=[]; fp(fp>1)=[];
%f=fp;%length(fp);
fpx=sortrows(fp);
fpz=arrayfun(@(X) xnull(X,a,g),fpx);
%}
%Plot setup:
%u=1; v=1;
u=1; v=1/a;
eps=.01;
X=(0:eps:u);
Z=(0:eps/a:v);%check with k in contumlim & at bottom
lx=length(X);
lz=length(Z);
%Comment out if input A:
%
%Basins:
%pm=(a-sqrt(a^2-4*a*v))/(2*a); pp=(a+sqrt(a^2-4*a*v))/(2*a);
%qm=(b-sqrt(b^2-4*b*k*u))/(2*b*k); qp=(b+sqrt(b^2-4*b*k*u))/(2*b*k);
xoverFcn=@(t,Y)evUN(t,Y,u,v); 
options=odeset('Events',xoverFcn);
if basins==1
    A=zeros(lx,lz);
    epsilon=.001;
    for i=1:lx
        parfor j=1:lz
            [t,y]=ode23(@(t,Y)integr8(t,Y,a,g,b1,b2,u,v),[0,2300],[X(i),Z(j)],options);
            %Find closest f.p.:
            
            fpdiff=sqrt((y(end,1)-fpx).^2+(y(end,2)-fpz).^2);
            [~,nfp]=min(fpdiff);
            A(i,j)=a*fpz(nfp)-fpx(nfp);
            
            %{
            if u-y(end,1)<epsilon %&& y(end,2)<epsilon
                A(i,j)=1;%cx light
            elseif v-y(end,2)<epsilon/a %y(end,1)<epsilon %&& v-y(end,2)<epsilon
                A(i,j)=2;%cz light
            end
            %}
        end
    end
end
f=A;
%}
%Figure:
fs=20; lw=2; ms=10;
figure
cx=[0,0,.5]; cz=[.5,0,0];
%map=[1,1,1;0.5,.5,1;0.5,.5,1;1,.5,.5];
%colormap(map)
colormap(.5+.5*redblue)
if basins>0
    %contourf(X,Z,A',3,'LineStyle','none');%'LineWidth',2);%,'LineStyle','none');%,'LineWidth',2)
    imagesc(X,Z,A')
    caxis([-1,1])
    set(gca,'YDir','normal')
end
hold on
xx=(0:.005:1);
zz=(0:.005*1/a:1/a);
p1=arrayfun(@(X) xnull(X,a,g),xx);
p2=arrayfun(@(Z) znull(Z,a,b1,b2),zz);
q1=arrayfun(@(X) Xquad(X,a),xx);
q2=arrayfun(@(Z) Yquad(Z,a,b1),zz);
plot(xx,p1,'color',cx,'LineWidth',lw)
%1 area in top right:
%plot(xx,q1,'--','color',cx,'LineWidth',lw)
if u<1
    plot([u,u],[0,1/k],'k--','linewidth',lw)
end
if v<1/a
    plot([0,1],[v,v],'k--','linewidth',lw)
end
plot (p2,zz,'color',cz,'LineWidth',lw)
%1 area in top right:
%plot (q2,zz,'--','color',cz,'LineWidth',lw)
%
%plot(u,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(0,v,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(.5,.5*a,'ks','markersize',ms,'markerfacecolor','k','linewidth',1)
%
%Fixed points (roots of p9 in order):
%Stable:
plot(fpx,fpz,'ko','markersize',ms,'markerfacecolor','k','linewidth',lw)
%Saddle:
sad=[2,4,6];
plot(fpx(sad),fpz(sad),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw)
plot(fpx(sad),fpz(sad),'kx','markersize',ms,'markerfacecolor','k','linewidth',lw)
%Unstable:
uns=5;
plot(fpx(uns),fpz(uns),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw)%[.75,.75,.75]

%Top right box:
%plot([.5,1],[1/2/a,1/2/a],'k-','LineWidth',.5)
%plot([.5,.5],[1/2/a,1/a],'k-','LineWidth',.5)
%
hold off
if u<1
    set(gca,'XTick',[0,u,1],'xticklabels',{'0','u','1'})
else
    set(gca,'XTick',[0,1])
end
if v<1/a
    set(gca,'YTick',[0,v,1/a],'yticklabels',{'0','v','1/{\alpha}'})
else
    set(gca,'YTick',[0,Z(end)],'yticklabels',{'0','1/{\alpha}'})
end
set(gca,'FontSize',fs)
xlabel('X_1')
ylabel('Y_1','rot',0)
axis([0,1,0,1/a])
caxis([-1,1])
%axis([.6,.75,.3/a,.35/a])
%grid on
%grid minor
box on
%{
%For case 1a: need b1=b2, gamma=1 here
figure
hold on
eps=.001;
xx=(1/2+eps:eps:1-eps);
zz=(1/2/a+eps/a:eps/a:1/a-eps/a);
%{
[XX,ZZ]=meshgrid(xx,zz);
[xx,zz]=XZtransform(xx,zz,a,b1);
[p1,p2]=XZtransform(p1,p2,a,b1);
[p1,p2]=XZtransform(xx,zz,a,b1);
%}
p1=arrayfun(@(X) xnull(X,a,g),xx);%So vectors have same length
p2=arrayfun(@(Z) znull(Z,a,b1,b2),zz);
q2=XZtransform(xx,p1,a,b1);
q1=XZtransform(p2,zz,a,b1);
plot(xx,q1,'color',cx,'LineWidth',lw)
plot (p2,q2,'color',cz,'LineWidth',lw)
hold off
set(gca,'FontSize',fs)
xlabel('X_1')
ylabel('Z_1','rot',0)
%axis([0,1,0,1/a])
%}
end

function f=xnull(X,a,g)
f=-(1+g)*X^3+(1+2*g)*X^2-(g+1/a)*X+1/a;
end

function g=znull(Z,a,b1,b2)
g=-a^2*(b1+b2)*Z^3+a*(b1+2*b2)*Z^2-(a+b2)*Z+1;
end

function [value,isterminal,direction]=evUN(t,Y,u,v)
l=length(Y);
value=Y-[u;v];
isterminal=ones(l,1);
direction=ones(l,1);
end

function f=integr8(t,Y,a,g,b1,b2,u,v)
%x=min(Y(1),u); y=min(Y(2),v);
X=Y(1); Z=Y(2);
f1=-(1+g)*X^3+(1+2*g)*X^2-(g+1/a)*X+1/a-Z;
f2=-a*(b1+b2)*Z^3+(b1+2*b2)*Z^2-(1+b2/a)*Z+1/a-X/a;
if X>=u && f1>0
    f1=0;
end
if Z>=v && f2>0
    f2=0;
end
f=[f1;f2];
end

function [f,g]=XZtransform(x,z,a,b)
f=1/a*(a*z+1-x)./x+1/a;
g=b*(x+1-a*z)./z+1/2;
%f=(1-x).*(2*x-1)+1/2/a;%2-2*x./(1-a*z)-1/2/a;
%g=.5*b*(1-a*z).*(2*z-1/a)+1/2;
%(z./(x-1)+1/a)*(x-1
end

function f=Xquad(x,a)
f=(2*x-1).*(1-x)+1/2/a;
end

function f=Yquad(z,a,b)
b=b*a/8;
f=b*(2*z-1/a).*(1-a*z)+1/2;
end