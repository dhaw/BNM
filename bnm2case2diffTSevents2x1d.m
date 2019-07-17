function f=bnm2case2diffTSevents2x1d(a,b1,x0,z0,plotTraj)%(a,b1,x0,z0,plotTraj)%(eps)
%Outputs 2 plots: 
%populations/time, and end state in phase space (with nullclines)
%Parameters:
%Same TS in both areas - gamma=1, beta2=beta1
%
%plotTraj=1;
%a=6; %alpha
g=1;
%b1=30; %beta1
b2=b1;
%}
tfinal=200;
%Initialconditions:
%{
x0=[.1,.11];
z0=[.09,.11]/a;
%}
%{
x0=rand(1,2); x0=x0/sum(x0);
z0=rand(1,2); z0=z0/sum(z0)/a;
%}
%%
y0=reshape([x0;z0],4,1);
tstart=0;
tout=tstart;
yout=y0';
teout=[];
yeout=[];

ieout=[];

options=odeset('Refine',1);
%options=odeset(options,'Events',@evUN);
%Integrate
[tout,yout]=ode23(@(t,Y)fsolve(t,Y,a,g,b1,b2),[tstart tfinal],y0,options);%,te,ye,ie
%}
%After integration
%Create output vectors for x1,z1,x2,z2
x1=yout(:,1); x2=yout(:,3);
z1=yout(:,2); z2=yout(:,4);
f=yout(end,:);
xneg=[x1,x2]; xneg(xneg>=0)=NaN;
zneg=[z1,z2]; zneg(zneg>=0)=NaN;

if plotTraj==1
%Plot
fs=15; lw=2; ms=10;
%cx=1.5*[0,.2,.2]; cz=1.5*[.149,.012,.224];
cx=[0,0,.5]; cz=[.5,0,0];
minxz=min(min(yout));
figure
subplot(2,1,1)%('position',[.05,.6,.45,.2])
hold on
h2=plot(tout,x1,'-','linewidth',lw,'color',cx);
h3=plot(tout,x2,'--','linewidth',lw,'color',cx);
plot(tout(end),x1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
plot(tout(end),x2(end),'kx','markersize',ms,'linewidth',lw);
plot(tout,xneg,'--m','linewidth',lw);
hold off
axis ([0,tfinal,minxz,1])
grid on
grid minor
box on
%xlabel('time','fontsize',15);
ylabel('X','rot',0)
set(gca,'fontsize',fs)
set(gca,'YTick',[0,1],'yticklabels',{'0','1'})
legend([h2,h3],'X_1','X_2','location','northeastoutside')
%
subplot(2,1,2)%('position',[.05,.2,.45,.2])
hold on
hh2=plot(tout,z1,'-','linewidth',lw,'color',cz);
hh3=plot(tout,z2,'--','linewidth',lw,'color',cz);
plot(tout(end),z1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
plot(tout(end),z2(end),'kx','markersize',ms,'linewidth',lw);
plot(tout,zneg,'--m','linewidth',lw);
hold off
axis ([0,tfinal,minxz,1/a])
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('Y','rot',0)
set(gca,'fontsize',fs)
set(gca,'YTick',[0,1/a],'yticklabels',{'0','1/\alpha'})
legend([hh2,hh3],'Y_1','Y_2','location','northeastoutside')
hold off
end
%Option to comment out: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%If stable fixed point, plot in phase portrait
%(one dot for each variable - same set of nullclines):
xx=(0:.005:1); p1=xx.*(1-xx); q1=g*xx.*(1-xx);
zz=(0:.005*1/a:1/a); p2=b1*zz.*(1-a*zz); q2=b2*zz.*(1-a*zz);
%fs=20;
%
figure
hold on
h1=plot(xx,q1,'color',cx,'LineWidth',lw);
h2=plot (q2,zz,'color',cz,'LineWidth',lw);
h3=plot(x1(end),z1(end),'ko','markersize',ms,'markerfacecolor','k','linewidth',lw);
h4=plot(x2(end),z2(end),'kx','markersize',ms,'linewidth',lw);
hold off
set(gca,'FontSize',fs)
set(gca,'XTick',[0,1],'yticklabels',{'0','1'})
set(gca,'YTick',[0,1/a],'yticklabels',{'0','1/\alpha'})
xlabel('X')
ylabel('Y','rot',0)
axis([0,1,0,1/a])
grid on
grid minor
box on
legend([h3,h4],'a_1','a_2','location','northeastoutside')
%}
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
[xdot,zdot]=derivs(x,z,a,g,b1,b2);
f2=reshape([xdot';zdot'],4,1);
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