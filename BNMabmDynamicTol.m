function [f,g]=BNMabmDynamicTol
%Fixed tolerance
a=2; b=2; tend=200;
%Initial condition:
%Number in each cell
x0=[300,100];%;3*ones(100,1)];
y0=[100,300];%;3*ones(100,1)];
totx=sum(x0);
toty=sum(x0);
%
n=length(x0);

%Linear tolerance schedules:
%{
xtol=fliplr((0:a/(totx-1):a)); %reptolx=repmat(Xtol',1,2);
ytol=fliplr((0:b/(toty-1):b)); %reptoly=repmat(Ytol',1,2);
%}
%Not including zero:
%
xtol=fliplr((a/totx:a/totx:a)); %reptolx=repmat(Xtol',1,2);
ytol=fliplr((b/toty:b/toty:b)); %reptoly=repmat(Ytol',1,2);
%}
%Lunch time

X=zeros(tend+1,n); Y=X;
X(1,:)=x0; Y(1,:)=y0;
Xdiff=zeros(n,2); Ydiff=Xdiff;%Number that relocate at each time step
xx=x0; yy=y0;
for t=2:tend+1
    %X-population:
    threshx=xx./yy;%Agents included in own nbhd
    threshx(yy==0)=0;
    xrem=zeros(1,n);
    for i=1:n
        vxi=xtol(1:xx(i));
        xrem(i)=length(vxi(vxi<threshx(i)));
    end
    Xdiff(t-1,:)=xrem;
    xx=xx-xrem+fliplr(xrem);
    X(t,:)=xx;
    %Y-population:
    threshy=yy./xx;%Agents included in own nbhd
    threshy(xx==0)=0;
    yrem=zeros(1,n);
    for i=1:n
        vyi=ytol(1:yy(i));
        yrem(i)=length(vyi(vyi<threshy(i)));
    end
    Ydiff(t-1,:)=yrem;
    yy=yy-yrem+fliplr(yrem);
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
figure
subplot(2,1,1)
%subplot(2,2,1)
%colormap(cx)
hold on
h1=plot(1:tend+1,X(:,1),'-','linewidth',lw,'color',cx1);
h2=plot(1:tend+1,X(:,2),':','linewidth',lw,'color',cx2);
hold off
axis([0,tend+1,0,totx])
set(gca,'YTick',[0,totx],'yticklabels',{'0','X_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('X','rot',0)
set(gca,'fontsize',fs)
legend([h1,h2],'X_1','X_2','location','northeastoutside')
%
subplot(2,1,2)
%subplot(2,2,3)
%colormap(cy)
hold on
h1=plot(1:tend+1,Y(:,1),'-','linewidth',lw,'color',cy1);
h2=plot(1:tend+1,Y(:,2),':','linewidth',lw,'color',cy2);
hold off
axis([0,tend+1,0,toty])
set(gca,'YTick',[0,toty],'yticklabels',{'0','Y_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('Y','rot',0)
set(gca,'fontsize',fs)
legend([h1,h2],'Y_1','Y_2','location','northeastoutside')
%{
subplot(2,2,2)
%colormap(cx)
hold on
h1=plot(1:tend,Xdiff(:,1),'-','linewidth',lw,'color',cx1);
h2=plot(1:tend,Xdiff(:,2),':','linewidth',lw,'color',cx2);
hold off
axis([0,tend,0,totx])
set(gca,'YTick',[0,totx],'yticklabels',{'0','X_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('\Delta X','rot',0)
set(gca,'fontsize',fs)
legend([h1,h2],'X_1','X_2','location','northeastoutside')
%
subplot(2,2,4)
%colormap(cy)
hold on
g1=plot(1:tend,Ydiff(:,1),'-','linewidth',lw,'color',cy1);
g2=plot(1:tend,Ydiff(:,2),':','linewidth',lw,'color',cy2);
hold off
axis([0,tend,0,toty])
set(gca,'YTick',[0,toty],'yticklabels',{'0','Y_{max}'})
grid on
grid minor
box on
xlabel('time','fontsize',15); ylabel('\Delta Y','rot',0)
set(gca,'fontsize',fs)
legend([g1,g2],'Y_1','Y_2','location','northeastoutside')
%}