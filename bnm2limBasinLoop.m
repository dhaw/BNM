function f=bnm2limBasinLoop
%Eg2:
a=7; b=49;
u=[1,.89,.75,.25,.11,.09];
lu=length(u);
basinCell=cell(1,lu);
x=NaN;
%Manually input f.p.s on interior
%(numbered by real root of p_9(X):
sad=[2,4,6,8;
    2,4,6,8;
    2,4,6,x;
    2,x,x,x;
    x,x,x,x;
    x,x,x,x];
stab=[1,3,7,9;
    1,3,x,x;
    1,3,x,x;
    1,3,x,x;
    1,x,x,x;
    1,x,x,x];
uns=[5;
    5;
    5;
    x;
    x;
    x];
ustab=[1,x;
    1,x;
    1,3;
    x,x;
    1,x;
    x,x];
usad=[x;
    x;
    2;
    x;
    2;
    x];
for i=1:lu
    stabi=stab(i,:); stabi(isnan(stabi)==1)=[];
    sadi=sad(i,:); sadi(isnan(sadi)==1)=[];
    unsi=uns(i,:); unsi(isnan(unsi)==1)=[];
    ustabi=ustab(i,:); ustabi(isnan(ustabi)==1)=[];
    usadi=usad(i,:); usadi(isnan(usadi)==1)=[];
    Ai=bnm2basinLim(a,b,u(i),stabi,sadi,unsi,ustabi,usadi);
    basinCell{i}=Ai;
end
f=basinCell;
%{
%Eg1:
a=6; b=30;
u=[1,.8,.71,.6,.4,.29];
lu=length(u);
basinCell=cell(1,lu);
x=NaN;
%Manually input f.p.s on interior
%(numbered by real root of p_9(X):
sad=[2,4;
    2,4;
    2,4;
    2,x;
    2,x;
    x,x];
stab=[1,5;
    1,x;
    1,x;
    1,x;
    1,x;
    1,x];
uns=[3;
    3;
    3;
    3;
    x;
    x];
ustab=[1,x;
    1,x;
    1,3;
    1,x;
    1,x;
    x,x];
usad=[x;
    x;
    2;
    2;
    x;
    x];
%Eg2:
a=7; b=49;
u=[1,.89,.75,.25,.11,.09];
lu=length(u);
basinCell=cell(1,lu);
x=NaN;
%Manually input f.p.s on interior
%(numbered by real root of p_9(X):
sad=[2,4;
    2,4;
    2,4;
    2,x;
    2,x;
    x,x];
stab=[1,6;
    1,x;
    1,x;
    1,x;
    1,x;
    1,x];
uns=[3;
    3;
    3;
    3;
    x;
    x];
%}