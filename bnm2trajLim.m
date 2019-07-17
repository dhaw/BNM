function [f,g]=bnm2trajLim%(A)%params/y must match contunlim & plotnullclines
%Case 1a, limit x only - upper or upper and lower
a=7; b=49;
u=.9;
ub=1-u;%u bottom - <.5, u>.5
v=.9/a;
vb=1/a-v;
y0=[.49,.51/a];
%
eps=.005;%0.001 for pretty plots
X=(ub:eps:u);
Y=[(vb:eps:v)];%/2),(1/a/2+.01:.01:1/a)]; %check with k in contumlim & at bottom
lx=length(X);
ly=length(Y);
A=zeros(lx,ly);
epsilon=.001;

%pm=(a-sqrt(a^2-4*a*v))/(2*a); pp=(a+sqrt(a^2-4*a*v))/(2*a);
%qm=(b-sqrt(b^2-4*b*k*u))/(2*b*k); qp=(b+sqrt(b^2-4*b*k*u))/(2*b*k);

xoverFcn=@(t,Y)evUN(t,Y,u,ub,v,vb); 
options=odeset('Events',xoverFcn);

tout=[];
yout=[];
teout=[];
yeout=[];
ieout=[];
tstart=0;
tend=200;
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
f=tout;
g=yout;
%%
fs=20; lw=2; ms=10;
xvec=yout(:,1); yvec=a*yout(:,2);
cx=[0,0,.5]; cy=[.5,0,0];
figure
hold on
h2=plot(tout,xvec,'-','linewidth',lw,'color',cx);
h3=plot(tout,yvec,'--','linewidth',lw,'color',cy);
%plot(tout(end),x1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
%plot(tout(end),x2(end),'kx','markersize',ms,'linewidth',lw);
%plot(tout,xneg,'--m','linewidth',lw);
hold off
axis ([0,tout(end),0,1])
grid on
grid minor
box on
xlabel('time','fontsize',12);
ylabel('population','fontsize',12);%,'rot',0
set(gca,'fontsize',fs)
set(gca,'YTick',[0,1],'yticklabels',{'0','1'})
legend([h2,h3],'X_1','aY_1','location','NE')
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
if y<=v && f2<0
    f2=0;
end
f=[f1;f2];
end