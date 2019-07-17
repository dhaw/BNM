function [f1,f2,numnans]=HPCabnUN
nn=2;%[2,4,6,10,20];
ln=length(nn);
a=2;%(2:2:20);%(1:50);
la=length(a);
%t=2*ones(la,1);
%t(1:5)=10*ones(5,1);
N=1;
numnans=0;

AA=zeros(la,ln,N);%1st arg=number of lines to plot
%BB=AA;
for I=1:N
for i=1:la%5 diff values of n/m
for j=1:ln
   % temp=zeros(N,1);
    for n=1:N
        try
        dissim=oneModel(a(i),a(i),nn(j),10);
        dissim=dissim(3);%Change for homoph x/y
        %if isnan(dissim)==0
            AA(i,j,I)=dissim;
            %BB(i,j)=BB(i,j)+1;
        %else
            %AA(i,n)=1;
        %end
        %AA(i,j,n)=dissim;
        %temp(n)=dissim;
        catch ME
            numnans=numnans+1;
            AA(i,j,I)=NaN;
        end
    end
    %AA(i,j)=AA(i,j)/BB(i,j);
   %A(i,j)=nanmean(temp);
end
end
end
f1=AA;
f2=nanmean(f1,3);
%save('abnUN.mat','f1','f2','numnans')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction]=evUN(t,y)
l=length(y);
value=y;
isterminal=ones(l,1);
direction=zeros(l,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f1=oneModel(a,b,n,tt)
k=4;
tfinal=tt;%Enough to guarantee an event?
%Adj mat:
A=ones(n,n);%Complete
A=A-eye(n);
%Random initial conditions:
x0=rand(1,n);
x0=x0*n/sum(x0);
y0=rand(1,n);
y0=y0*n/(k*sum(y0));
xy=[x0;y0];
y0=reshape(xy,1,2*n);
options=odeset('Refine',1);
options=odeset(options,'Events',@evUN);
tstart=0;
tout=tstart;
yout=y0;
teout=[];
yeout=[];
ieout=[];
while tstart<tfinal
%Solve until the first terminal event:
[t,y,te,ye,ie]=ode45(@(t,Y)ode(t,Y,a,b,A,k),[tstart tfinal],y0,options);
%Accumulate output:  
nt=length(t);
tout=[tout;t(2:nt)];
yout=[yout;y(2:nt,:)];
teout=[teout;te];
yeout=[yeout;ye];
ieout=[ieout;ie]; 
%Set the new initial conditions:
y0=y(nt,:);
y0(y0<=0)=0;
    xx=y0(1:2:end-1);
    yy=y0(2:2:end);
    for i=1:length(ie);
    if mod(ie(i),2)==0
    	index=ie(i)/2;
    	dy=yy(1,index);
    	[cy,Iy]=max(yy);
    	Iy=Iy(1);
    	yy(Iy)=yy(Iy)+dy;
     	yy(index)=0;
    else
    	index=ceil(ie(i)/2);
    	dx=xx(1,index);
    	[cx,Ix]=max(xx);
    	Ix=Ix(1);
     	xx(Ix)=xx(Ix)+dx;
    	xx(index)=0;
    end
    end
xy=[xx;yy];
y0=reshape(xy,1,2*n);   
yout(end,:)=y0;%Puts new zero in right place
tstart=t(nt);
end
%X=zeros(n,2);
%Checking conservation:
x=yout;
xx=0;
yy=0;
g=zeros(n,2);
for i=1:n
    g(i,1)=x(end,2*i-1);
    g(i,2)=x(end,2*i);
    xx=sum(g(:,1));
    yy=sum(g(:,2));
end
gg=[xx,yy];
gx=x(end,:);
%for i=1:n
%    X(i,1)=gx(2*i-1);
%    X(i,2)=gx(2*i);
%end
X=reshape(gx,2,n);
X=X';
xx=X(:,1);
yy=X(:,2);
nn=xx+yy;
%d=zeros(n,1);
%HX=d;
%HY=d;
wx=k/(k+1);
wy=wx/k;
%{
for i=1:n
    d(i)=abs(g(i,1)-k*g(i,2));
    HX(i)=g(i,1)/(g(i,1)+g(i,2));
    HIX(i)=(HX(i)-wx)/(1-wx);
    HY(i)=g(i,2)/(g(i,1)+g(i,2));
    HIY(i)=(HY(i)-wy)/(1-wy);
end
%}
d=abs(g(i,1)-k*g(i,2));
dissim=.5/n*sum(d,2);
HX=xx./nn;
HY=yy./nn;
HIX=(HX-wx)/(1-wx);
HIY=(HY-wy)/(1-wy);

homophx=nanmean(HX);
homophy=nanmean(HY);
HIX(isnan(HIX))=0;%Indices change if entries removed - OK for computing MX as contributes 0
HIY(isnan(HIY))=0;
f1=[sum(X(:,1)),sum(X(:,2)),dissim,homophx,homophy];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f2=ode(t,Y,a,b,A,k)%Delta
n=length(Y)/2;
f2=zeros(2*n,1);
%F and G:
%Linear:
HF=@(xx,yy) -a*xx^3/n+a*xx^2-xx*yy;
HG=@(xx,yy) -k*b*yy^3/n+b*yy^2-xx*yy;
%Diff equations:
x=Y(1:2:end-1);
y=Y(2:2:end);
[BX,BY]=meshgrid(x,y);
BY=BY';
%BX=BX';
BX=BX.*A;
BY=BY.*A;
CX=BX./(BX+BY);
CY=BY./(BX+BY);
CX(isnan(CX))=0;
CY(isnan(CY))=0;
CX(CX>0)=1;%Comment out for PR
CY(CY>0)=1;%Comment out for PR
xfact=sum(CX,2);%Column vector-row sums
yfact=sum(CY,2);
xfact(xfact==0)=1;%Affects rest of calc?
yfact(yfact==0)=1;
P=zeros(n,n);
Q=P;
for i=1:n
        P(i,:)=CX(i,:)/xfact(i);
        Q(i,:)=CY(i,:)/yfact(i);
end
P(isnan(P))=0;
Q(isnan(Q))=0;
F=zeros(n,1);
G=F;
for i=1:n
	F(i)=HF(x(i),y(i));
	G(i)=HG(x(i),y(i));
end
for i=1:n
    if any(P(:,i))==0%nkx(whichcompx(i))==x(i)
        f2(2*i-1,1)=0;
    else
        f2(2*i-1,1)=F(i)-dot(P(:,i),F);
    end
    if any(Q(:,i))==0%nky(whichcompy(i))==y(i)
        f2(2*i,1)=0;
    else
        f2(2*i,1)=G(i)-dot(Q(:,i),G);
    end
end
end