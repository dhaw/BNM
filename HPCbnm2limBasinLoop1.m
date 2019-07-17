function f=HPCbnm2limBasinLoop1
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
for i=1:lu
    stabi=stab(i,:); stabi(isnan(stabi)==1)=[];
    sadi=sad(i,:); sadi(isnan(sadi)==1)=[];
    unsi=uns(i,:); unsi(isnan(unsi)==1)=[];
    ustabi=ustab(i,:); ustabi(isnan(ustabi)==1)=[];
    usadi=usad(i,:); usadi(isnan(usadi)==1)=[];
    Ai=HPCbnm2basinLim(a,b,u(i),stabi,sadi,unsi,ustabi,usadi);
    basinCell{i}=Ai;
end
f=basinCell;
save('basinLoopLim1.mat','basinCell')