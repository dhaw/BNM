function f=HPCbnm2basinLim(a,b,u,stab,ustab)
ub=0;%1-u;%u bottom - <.5, u>.5
v=1/a;
vb=1/a-v;
eps=.01;%For basins - 0.001 for pretty plots
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
fpx=sortrows(fp);
fpz=arrayfun(@(X) xnull(X,a,1),fpx);
%
%For upper limit u:
fpuz=roots([-2*a^2*b,3*a*b,-a-b,1-u]); fpuz=double(fpuz); fpuz=round(fpuz,4);
fpuz(imag(fpz)~=0)=[];
zupper=(1-u)*(1-a*u+2*a*u^2);
fpuz(fpuz<0)=[]; fpuz(fpuz>zupper)=[];
fpuz=sortrows(fpuz);
fpux=u*ones(length(fpuz),1);
stabfpx=[fpx(stab);fpux(ustab)];
stabfpz=[fpz(stab);fpuz(ustab)];
%}
%
X=(ub:eps:u);
Y=[(vb:eps/a:v)];%check with k in contumlim & at bottom
lx=length(X);
ly=length(Y);
A=zeros(lx,ly);
xoverFcn=@(t,Y)evUN(t,Y,u,ub,v,vb); 
options=odeset('Events',xoverFcn);

tout=[];
yout=[];
teout=[];
yeout=[];
ieout=[];
tend=1000;
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
        y0=y(end,:);
        tstart=t(end);
        end
        %
        %Find closest f.p.:
        fpdiff=sqrt((y(end,1)-stabfpx).^2+(y(end,2)-stabfpz).^2);
        [~,nfp]=min(fpdiff);
        aij=a*stabfpz(nfp)-stabfpx(nfp);
        yveci(j)=aij;%A(i,j)=aij;
    end
    A(i,:)=yveci;
end
f=A;
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

function f=xnull(X,a,g)
f=-(1+g)*X^3+(1+2*g)*X^2-(g+1/a)*X+1/a;
end