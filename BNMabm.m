function [f,g]=BNMabm
randTol=1;
%Fixed tolerance
%To do: include diff (number that relocate - may be more than simple difference!)
a=2; b=2; tend=200;
%Initial condition:
%1 means "in cell 1" etc.
Xcell=[ones(300,1);2*ones(100,1)];%;3*ones(100,1)];
Ycell=[ones(100,1);2*ones(300,1)];%;3*ones(100,1)];
%
lx=length(Xcell);
ly=length(Ycell);
n=max([Xcell;Ycell]);
k=lx/ly;
nnx=repmat(cumsum(ones(1,n)),lx,1); nny=repmat(cumsum(ones(1,n)),ly,1);
Nx=repmat(Xcell,1,n); Ny=repmat(Ycell,1,n);
Nx(Nx~=nnx)=0; Nx(Nx>0)=1; Ny(Ny~=nny)=0; Ny(Ny>0)=1;
%Linear tolerance schedules:
Xtol=(a/lx:a/lx:a)'; 
if randTol==1
    Xtol=Xtol(randperm(lx));
end
Xtol=repmat(Xtol,1,n); Xtol=Xtol.*Nx;%No zero - (0:a/(lx-1):a)'
Ytol=(b/ly:b/ly:b)'; 
if randTol==1
    Ytol=Ytol(randperm(ly));
end
Ytol=repmat(Ytol,1,n); Ytol=Ytol.*Ny;
%Objects:
%one column per nbhd,
%IF agent present, then tolerance
%ELSE zero

%Nx=sparse(Nx); Ny=sparse(Ny);
%Xtol=sparse(Xtol); Ytol=sparse(Ytol);
X=zeros(tend+1,n); Y=X;
x0=sum(Nx,1); y0=sum(Ny,1);
X(1,:)=x0; Y(1,:)=y0;
Xdiff=zeros(tend,n); Ydiff=Xdiff;
for t=2:tend+1
    hx=sum(Nx,1); hy=sum(Ny,1);
    %X-population:
    threshx=hx./hy; threshx(hy==0)=0;%Agents included in own nbhd
    [agentx,cx,tolx]=find(Xtol); X3=[agentx,cx,tolx]; X3=sortrows(X3,2);%Row/col/val
    vx=repelem(threshx',hx);%Vector - one entry per agent, cell frac, *ordered by cell*
    movex=tolx-vx; movex(movex>0)=0; movex(movex<0)=1;%Which agents to move
    movex=movex.*X3(:,2);%Where from
    findx=find(movex); lfx=length(findx); rx=randsample(n-1,lfx,true);%Move at random
        
    if isempty(findx)==0
        cxfind=cx(findx);
        Xdifft=[length(cxfind(cxfind==1)),length(cxfind(cxfind==2))];%accumarray(cx(findx),ones(lfx,1));
        Xdiff(t-1,:)=Xdifft;
        %
        cxto=mod(cx(findx)+rx,n);%Where to move to 
        cx(findx)=cxto; cx(cx==0)=n;%Re-index columns
        Xtol=sparse(agentx,cx,tolx,lx,n);%New loc/tol matrix
        Nx=Xtol; Nx(Nx>0)=1;
    end
    xx=sum(Nx,1);
    X(t,:)=xx;
    %Y-population:
    threshy=hy./hx; threshy(hx==0)=0;
    [agenty,cy,toly]=find(Ytol); Y3=[agenty,cy,toly]; Y3=sortrows(Y3,2);
    vy=repelem(threshy',hy);
    movey=toly-vy; movey(movey>0)=0; movey(movey<0)=1; 
    movey=movey.*Y3(:,2);
    findy=find(movey); lfy=length(findy); ry=randsample(n-1,lfy,true);
    
    if isempty(findy)==0
        cyfind=cy(findy);
        Ydifft=[length(cyfind(cyfind==1)),length(cyfind(cyfind==2))];%accumarray(cy(findy),ones(lfy,1));
        Ydiff(t-1,:)=Ydifft;
        %
        cyto=mod(cy(findy)+ry,n); 
        cy(findy)=cyto; cy(cy==0)=n;
        Ytol=sparse(agenty,cy,toly,ly,n); 
        Ny=Ytol; Ny(Ny>0)=1;
    end
    yy=sum(Ny,1);
    Y(t,:)=yy;
    %{
    bar([xx',yy']);%,'color',[0,0,.5]);%,'linewidth',2) ,'histc'
    axis([.5,n+.5,0,lx])
    xlabel('Cell','fontsize',15); ylabel('Population','fontsize',10)
    legend('X','Y')
    pause(0.5)%drawnow
    %}
end
f=X; g=Y;
%
%Plot
fs=15; lw=2; ms=10;
%cx=1.5*[0,.2,.2]; cz=1.5*[.149,.012,.224];
cx1=[0,0,.5]; cy1=[.5,0,0];
cx2=[0,0,.8]; cy2=[.8,0,0];
%scaleCol=repmat((1/n:1/n:1)',1,3);
%cx=repmat(cx,n,1).*scaleCol;
%cy=repmat(cy,n,1).*scaleCol;
%FigH=figure('DefaultAxesPosition', [0.1, 0.1, 0.9, 0.9]);
figure
subplot(2,1,1)
%subplot(2,2,1)
hold on
h1=plot(1:tend+1,X(:,1),'-','linewidth',lw,'color',cx1);
h2=plot(1:tend+1,X(:,2),':','linewidth',lw,'color',cx2);
hold off
axis([0,tend+1,0,lx])
set(gca,'YTick',[0,lx],'yticklabels',{'0','X_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('X','rot',0)
set(gca,'fontsize',fs)
legend([h1,h2],'X_1','X_2','location','northeastoutside')
%
subplot(2,1,2)
%subplot(2,2,3)
hold on
g1=plot(1:tend+1,Y(:,1),'-','linewidth',lw,'color',cy1);
g2=plot(1:tend+1,Y(:,2),':','linewidth',lw,'color',cy2);
hold off
axis([0,tend+1,0,ly])
set(gca,'YTick',[0,ly],'yticklabels',{'0','Y_{max}'})
grid on
grid minor
box on;
xlabel('time','fontsize',15); ylabel('Y','rot',0)
set(gca,'fontsize',fs)
legend([g1,g2],'Y_1','Y_2','location','northeastoutside')
%{
subplot(2,2,2)
hold on
h1=plot(1:tend,Xdiff(:,1),'-','linewidth',lw,'color',cx1);
h2=plot(1:tend,Xdiff(:,2),':','linewidth',lw,'color',cx2);
hold off
axis([0,tend,0,lx])
set(gca,'YTick',[0,lx],'yticklabels',{'0','X_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('\Delta X','rot',0)
set(gca,'fontsize',fs)
legend([h1,h2],'X_1','X_2','location','northeastoutside')
%
subplot(2,2,4)
hold on
g1=plot(1:tend,Ydiff(:,1),'-','linewidth',lw,'color',cy1);
g2=plot(1:tend,Ydiff(:,2),':','linewidth',lw,'color',cy2);
hold off
axis([0,tend,0,ly])
set(gca,'YTick',[0,ly],'yticklabels',{'0','Y_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('\Delta Y','rot',0)
set(gca,'fontsize',fs)
legend([g1,g2],'Y_1','Y_2','location','northeastoutside')
%}