function [f,g]=bnm2case2lyap(a,b1,x0,z0)
%Outputs 2 plots: 
%populations/time, and end state in phase space (with nullclines)
%Parameters:
%Same TS in both areas - gamma=1, beta2=beta1
%
%a=4; %alpha
g=1; %gamma=a2/a1
%b1=60; %beta1
b2=30; %beta2
zeros2=zeros(2);
%}
tfinal=100;
ltimestep=1;
%Initialconditions:
%{
x0=[.1,.3,.6]; x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
z0=[.4,.2,.4]; z0=z0/sum(z0)/a; %[z1,z2,z3] - can input as ratio
%}
%{
x0=rand(1,3); x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
z0=rand(1,3); z0=z0/sum(z0)/a; %[z1,z2,z3] - can input as ratio
%}
%%
y0=reshape([x0(1:2);z0(1:2)],4,1);
tstart=0;

options=odeset('Refine',1);
%options=odeset(options,'Events',@evUN);

[T,Res]=lyapunov2(4,@(t,Y)flyap(t,Y,a,g,b1,b2,zeros2),@ode45,0,ltimestep,tfinal,y0',0); 

f=T;
g=Res;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=flyap(t,Y,a,g,b1,b2,zeros2)
f=zeros(20,1);
f(1:4)=fsolve(t,Y(1:4),a,g,b1,b2);
YY=reshape(Y(5:end),4,4);
J=Jac(Y,a,g,b1,b2,zeros2);
f(5:end)=J*YY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f2=fsolve(t,Y,a,g,b1,b2)
x=Y([1,3]);
z=Y([2,4]);
[xdot,zdot]=derivs(x,z,a,g,b1,b2);
if 1-sum(x)<=0
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
else
    dx=xdot;
end
if 1/a-sum(z)<=0
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
else
    dz=zdot;
end

f2=reshape([dx';dz'],4,1);
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

function f=Jac(Y,a,g,b1,b2,zeros2)
X1=Y(1); X2=Y(3); Z1=Y(2); Z2=Y(4);
Phik1=1-sig(X1^2-X1^3-X1*Z1)*sig(X1+X2);
Phik2=1-sig(b1*Z1^2-a*b1*Z1^3-X1*Z1)*sig(a*Z1+a*Z2);
Phil1=1-sig(g*X2^2-g*X2^3-X2*Z2)*sig(X1+X2);
Phil2=1-sig(b2*Z2^2-a*b2*Z2^3-X2*Z2)*sig(a*Z1+a*Z2);
dk11=(2*X1-3*X1^2-Z1)*Phik1;
dk12=-X1*Phik1;
dk21=-Z1*Phik2;
dk22=(2*b1*Z1-a*b1*Z1^2-X1)*Phik2;
dl11=(2*g*X2-3*g*X2^2-Z2)*Phil1;
dl12=-X2*Phil1;
dl21=-Z2*Phil2;
dl22=(2*b2*Z2-a*b2*Z2^2-X2)*Phil2;
K=[dk11,dk12;dk21,dk22];
L=[dl11,dl12;dl21,dl22];
f=[K,zeros2;zeros2,L];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f=sig(X)
if X<0
    f=0;
else
    f=1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction]=evUN(t,y)
value=y(end-1:end);
isterminal=zeros(2,1);
direction=-ones(2,1);
end