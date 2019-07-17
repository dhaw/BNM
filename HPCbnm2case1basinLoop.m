function f=HPCbnm2case1basinLoop
%Eg1:
a=9; u=1;
b=[16,40,80];
lb=length(b);
basinCell=cell(1,lu);
x=NaN;
%Manually input f.p.s on interior
%(numbered by real root of p_9(X):
sad=[2,4;
    2,4;
    x,x];
stab=[1,5;
    1,x;
    1,x];
uns=[3;
    3;
    x];
for i=1:lb
    stabi=stab(i,:); stabi(isnan(stabi)==1)=[];
    sadi=sad(i,:); sadi(isnan(sadi)==1)=[];
    unsi=uns(i,:); unsi(isnan(unsi)==1)=[];
    Ai=HPCbnm2basinLim(a,b(i),u,stabi,sadi,unsi,[],[]);
    basinCell{i}=Ai;
end
f=basinCell;
save('basinLoopCase1.mat','basinCell')