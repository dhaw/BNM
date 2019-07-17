function f=bnm2case2diffTS(a,b1,x0,z0,plotTraj)%(eps)
%Outputs 2 plots: 
%populations/time, and end state in phase space (with nullclines)
%Parameters:
%Same TS in both areas - gamma=1, beta2=beta1
%{
a=3.1; %alpha
g=1; %gamma=a2/a1
b1=9.2; %beta1
b2=9.6; %beta2
%}
%
%a=6.1; %alpha
g=1;%4; %gamma=a2/a1
%b1=3; %beta1
b2=b1;%40; %beta2
%}
%{
a=3+3*rand; %alpha
g=0.5+3.5*rand; %gamma=a2/a1
b1=9+6*rand; %beta1
b2=9+6*rand; %beta2
%}
tfinal=100;
%Initialconditions:
%{
x0=[.1,.3,.6]; x0=x0/sum(x0); %[x1,x2,x3], normalised just in case [.1,.3,.6]
z0=[.4,.2,.4]; z0=z0/sum(z0)/a; %[z1,z2,z3] - can input as ratio
%}
%{
x0=[00.0117,0.5629,0.4254]; x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
z0=[0.0907,0.1095,0.1223]; z0=z0/sum(z0)/a; %[z1,z2,z3] - can input as ratio
%}
%{
x0=rand(1,3); x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
z0=rand(1,3); z0=z0/sum(z0)/a; %[z1,z2,z3] - can input as ratio
%}
%%
y0=reshape([x0;z0],6,1);
tstart=0;
tout=tstart;
yout=y0';
teout=[];
yeout=[];

ieout=[];

options=odeset('Refine',1);
options=odeset(options,'Events',@evUN);
%Integrate

while tstart<tfinal
%Solve until the first terminal event:
[t,y,te,ye,ie]=ode45(@(t,Y)fsolve(t,Y,a,g,b1,b2),[tstart tfinal],y0,options);%,te,ye,ie
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
y0(y0<0)=0;
y0=reshape(y0,2,3);
y01=y0(1,:); y02=y0(2,:);
%R-allocate according to existing proportions:
y01=y01/sum(y01);
y02=y02/sum(y02)/a;
y0=[y01;y02];
y0=reshape(y0,6,1);
tstart=t(end);
end
%}
%After integration
%Create output vectors for x1,z1,x2,z2
x1=yout(:,1); x2=yout(:,3); x3=yout(:,5);
z1=yout(:,2); z2=yout(:,4); z3=yout(:,6);
f=yout;
%f=[x0';z0';a;g;b1;b2];
%x3ind=find(ieout==5); lx3=length(x3ind);
%z3ind=find(ieout==6); lz3=length(z3ind);
xneg=[x1,x2,x3]; xneg(xneg>=0)=NaN;
zneg=[z1,z2,z3]; zneg(zneg>=0)=NaN;

if plotTraj==1
%Plot
fs=15; lw=2; ms=10;
%cx=1.5*[0,.2,.2]; cz=1.5*[.149,.012,.224];
cx=[0,0,.5]; cz=[.5,0,0];
minxz=min(min(yout));
figure
subplot(2,1,1)%('position',[.05,.6,.45,.2])
hold on
h1=plot(tout,x3,'-','linewidth',lw,'color','k');
h2=plot(tout,x1,'-','linewidth',lw,'color',cx);
h3=plot(tout,x2,'--','linewidth',lw,'color',cx);
plot(tout(end),x1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
plot(tout(end),x2(end),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw);
plot(tout(end),x3(end),'kx','markersize',ms,'linewidth',lw);
plot(tout,xneg,'--m','linewidth',lw);
%{
if lx3>0
    plot(teout(x3ind),[zeros(lx3,1),ones(lz3,1)],'--','k','linewidth',lw)
end
%}
hold off
axis ([0,tfinal,minxz,1])
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('X','rot',0)
set(gca,'fontsize',fs)
set(gca,'YTick',[0,1],'yticklabels',{'0','1'})
legend([h2,h3,h1],'X_1','X_2','X_3','location','northeastoutside')
%
subplot(2,1,2)%('position',[.05,.2,.45,.2])
hold on
hh1=plot(tout,z3,'-','linewidth',lw,'color','k');
hh2=plot(tout,z1,'-','linewidth',lw,'color',cz);
hh3=plot(tout,z2,'--','linewidth',lw,'color',cz);
plot(tout(end),z1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
plot(tout(end),z2(end),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw);
plot(tout(end),z3(end),'kx','markersize',ms,'linewidth',lw);
plot(tout,zneg,'--m','linewidth',lw);
hold off
axis ([0,tfinal,minxz,1/a])
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('Z','rot',0)
set(gca,'fontsize',fs)
set(gca,'YTick',[0,1/a],'yticklabels',{'0','1/\alpha'})
legend([hh2,hh3,hh1],'Z_1','Z_2','Z_3','location','northeastoutside')
hold off
end
%%
%Option to comment out: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%If stable fixed point, plot in phase portrait
%(one dot for each variable - same set of nullclines):
xx=(0:.005:1); p1=xx.*(1-xx); q1=g*xx.*(1-xx);
zz=(0:.005*1/a:1/a); p2=b1*zz.*(1-a*zz); q2=b2*zz.*(1-a*zz);
%fs=20;
%{
figure%subplot('position',[.55,.3,.35,.45])%figure
hold on
h1=plot(xx,p1,'color',cx,'LineWidth',lw);
h2=plot (p2,zz,'color',cz,'LineWidth',lw);
h1a=plot(xx,q1,'--','color',cx,'LineWidth',lw);
h2a=plot (q2,zz,'--','color',cz,'LineWidth',lw);
h3=plot(x1(end),z1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
h4=plot(x2(end),z2(end),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw);
h5=plot(x3(end),z3(end),'kx','markersize',ms,'linewidth',lw);
hold off
set(gca,'FontSize',fs)
set(gca,'XTick',[0,1],'yticklabels',{'0','1'})
set(gca,'YTick',[0,1/a],'yticklabels',{'0','1/\alpha'})
xlabel('X')
ylabel('Z','rot',0)
axis([minxz,1,minxz,1/a])
grid on
grid minor
box on
legend([h3,h4,h5],'a_1','a_2','a_3','location','northeastoutside')
%}
%{
figure
hold on
h1=plot(xx,q1,'color',cx,'LineWidth',lw);
h2=plot (q2,zz,'color',cz,'LineWidth',lw);
h3=plot(x1(end),z1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
h4=plot(x2(end),z2(end),'ko','markersize',ms,'markerfacecolor','w','linewidth',lw);
h5=plot(x3(end),z3(end),'kx','markersize',ms,'linewidth',lw);
hold off
set(gca,'FontSize',fs)
set(gca,'XTick',[0,1],'yticklabels',{'0','1'})
set(gca,'YTick',[0,1/a],'yticklabels',{'0','1/\alpha'})
xlabel('X')
ylabel('Z','rot',0)
axis([minxz,1,minxz,1/a])
grid on
grid minor
box on
legend([h3,h4,h5],'a_1','a_2','a_3','location','northeastoutside')
%}
%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction]=evUN(t,y)
value=y(end-1:end);
isterminal=zeros(2,1);%ones zeros
direction=-ones(2,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fx,fz]=derivs(x,z,a,g,b1,b2)
x1=x(1); x2=x(2);
z1=z(1); z2=z(2);
fx1=x1.^2.*(1-x1)-x1.*z1;
fx2=g*x2.^2.*(1-x2)-x2.*z2;
fz1=b1*z1.^2.*(1-a*z1)-x1.*z1;
fz2=b2*z2.^2.*(1-a*z2)-x2.*z2;
fx=[fx1;fx2]; fz=[fz1;fz2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f2=fsolve(t,Y,a,g,b1,b2)
x=Y([1,3]);
z=Y([2,4]);
xz3=Y([5,6]);
[xdot,zdot]=derivs(x,z,a,g,b1,b2);
if xz3(1)<=0
    if xdot(1)<0
        dx1=xdot(1);
    else
        dx1=0;
    end
    if xdot(2)<0
        dx2=xdot(2);
    else
        dx2=0;
    end
    dx=[dx1;dx2];
    %dx3=-sum(dx);;
else
    dx=xdot;
    %dx3=-sum(dx);
end
if xz3(2)<=0
    if zdot(1)<0
        dz1=zdot(1);
    else
        dz1=0;
    end
    if zdot(2)<0
        dz2=zdot(2);
    else
        dz2=0;
    end
    dz=[dz1;dz2];
    %dz3=-sum(dz);
else
    dz=zdot;
    %dz3=-sum(dz);
end
dx3=-sum(dx);
dz3=-sum(dz);
f2=reshape([dx',dx3;dz',dz3],6,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Stop on event - set the new initial conditions:
%{
y0=y(nt,:);
y0(y0<=0)=0;
%Just in case, one way to re-normalise: ********
xall=y0([1,3,5]);
%xall(xall>1)=1;
xall=xall/sum(xall);
zall=y0([2,4,6]);
%zall(zall>1/a)=1/a;
zall=zall/sum(zall);
xz=[xall;zall];
y0=reshape(xz,6,1);   
yout(end,:)=y0';%Puts new zero in right place
tstart=t(nt);
end
%}