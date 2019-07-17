function f=Zbasin%(A)%params/y must match contunlim & plotnullclines
a=1; b=9.2; k=3.1;
%u=1; v=.5;
u=1; v=1/k;
eps=.01;%0.001 for pretty plots
X=(0:eps:u);
Y=(0:eps:min(1/k,v)); %check with k in contumlim & at bottom
lx=length(X);
ly=length(Y);
A=zeros(lx,ly);
epsilon=.05;

pm=(a-sqrt(a^2-4*a*v))/(2*a); pp=(a+sqrt(a^2-4*a*v))/(2*a);
qm=(b-sqrt(b^2-4*b*k*u))/(2*b*k); qp=(b+sqrt(b^2-4*b*k*u))/(2*b*k);

xoverFcn=@(t,Y)evUN(t,Y,u,v); 
options=odeset('Events',xoverFcn);

for i=1:lx
    for j=1:ly
        [t,x]=ode23(@(t,Y)Zcontunlim(t,Y,a,b,k,u,v),[0,200],[X(i),Y(j)],options);
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
    end
end
f=A;

ali=a*k; betj=a*b;
xij=roots([1,-2,(1+ali)/ali,(1-betj)/(ali*betj)]); %xij(imag(xij)~=0)=[];
xij=sortrows(xij); yij=xij.*(1-xij)*a;
if length(xij)==3
    xo=xij(1); xo2=xij(3); xf=xij(2);
    yo=yij(1); yo2=yij(3); yf=yij(2);
else
    xo=xij(1); yo=yij(1);
end

%}
fs=20; lw=2; ms=10;
figure
map=[1,.5,.5;1,1,1;1,1,1;1,1,1;.5,.5,1];
%map=[1,.5,.5;1,1,1;1,1,1;1,1,1;.5,.5,1;.9,.9,.9];
%%map=[0,0,0;1,.5,.5;1,1,1;.5,.5,1;.9,.9,.9;1,1,1];
colormap(map)
contourf(X,Y,A',5,'LineStyle','none');%'LineWidth',2);%,'LineStyle','none');%,'LineWidth',2)%
axis([0,1,0,1/k])
hold on
xx=(0:.005:1);
yy=(0:.005:1/k);
[p1,p2]=Zplotnullclines(a,b,k,xx,yy);
plot(xx,p1,'color',[0 0 .5],'LineWidth',lw)
if u<1
    plot([u,u],[0,1/k],'k--','linewidth',lw)
end
if v<1/k
    plot([0,1],[v,v],'k--','linewidth',lw)
end
plot (p2,yy,'color',[.5 0 0],'LineWidth',lw)
plot([0,1],[0,0],'color',[.5 0 0],'LineWidth',lw)
plot([0,0],[0,1/k],'color',[0 0 .5],'LineWidth',lw)

%
plot(u,0,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
plot(0,v,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)

plot(xo,yo,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(xf,yf,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(xo2,yo2,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(pm,v,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)
%plot(u,qm,'ko','markersize',ms,'markerfacecolor','w','linewidth',1)

%plot(u,a*u*(1-u),'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%plot(b*v*(1-k*v),v,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)

%
hold off
if u<1
    set(gca,'XTick',[0,u,1],'xticklabels',{'0','u','1'})
else
    set(gca,'XTick',[0,1])
end
if v<1/k
    set(gca,'YTick',[0,v,1/k],'yticklabels',{'0','v','0.5'})
else
    set(gca,'YTick',[0,Y(end)],'yticklabels',{'0','0.5'})
end
set(gca,'FontSize',fs)
xlabel('X')
ylabel('Y','rot',0)
axis([0,1,0,1/k])
%f=x;
end

function [value,isterminal,direction]=evUN(t,Y,u,v)
l=length(Y);
value=Y-[u;v];
isterminal=ones(l,1);
direction=ones(l,1);
end

function f=Zcontunlim(t,Y,a,b,k,u,v)
%x=min(Y(1),u); y=min(Y(2),v);
x=Y(1); y=Y(2);
f1=a*x^2*(1-x)-y*x;
f2=b*y^2*(1-k*y)-x*y;
if x>=u && f1>0
    f1=0;
end
if y>=v && f2>0
    f2=0;
end
f=[f1;f2];
end

function [f,g]=Zplotnullclines(a,b,k,X,Y)
f=a*X.*(1-X);
g=b*Y.*(1-k*Y);
end