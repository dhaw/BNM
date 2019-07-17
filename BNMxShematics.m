function f=BNMxShematics
eps=.05;
tols=fliplr(0:eps:1);
lx=length(tols);
xvec=(1:lx);

fs=12;
lw=2;
ms=10;
%{
figure
bar(xvec,tols)
axis([.5,xvec(end)+.5,0,1])
set(gca,'XTick',[xvec(1),xvec(end)],'xticklabels',{'0','1'})
set(gca,'YTick',[0,1],'yticklabels',{'0','a'})
grid on
grid minor
box on
xlabel('X','fontsize',15); ylabel('R_X(X)')%,'rot',0)
set(gca,'fontsize',fs)
%legend([h1,h2],'X_1','X_2','location','northeastoutside')
%}
%
thismany=round(.4*lx);
x1=randsample(lx,thismany);
bar1=tols; bar1(x1)=0;
bar2=tols-bar1;

figure
subplot(2,2,1)
hold on
bar(xvec(1:thismany),tols(1:thismany))
bar(xvec(thismany+1:end),tols(thismany+1:end),'facecolor',.7*[1,1,1])
hold off
axis([.5,xvec(end)+.5,0,1])
set(gca,'XTick',[xvec(1),xvec(end)],'xticklabels',{'0','1'})
set(gca,'YTick',[0,1],'yticklabels',{'0','a'})
grid on
grid minor
box on
xlabel('X_1','fontsize',15); ylabel('R_X(X_1)')%,'rot',0)
set(gca,'fontsize',fs)
title('Nbhd 1')
%legend([h1,h2],'X_1','X_2','location','northeastoutside')

subplot(2,2,2)
hold on
bar(xvec(1:lx-thismany+1),tols(1:lx-thismany+1))
bar(xvec(lx-thismany+1:end),tols(lx-thismany+1:end),'facecolor',.7*[1,1,1])
hold off
axis([.5,xvec(end)+.5,0,1])
set(gca,'XTick',[xvec(1),xvec(end)],'xticklabels',{'0','1'})
set(gca,'YTick',[0,1],'yticklabels',{'0','a'})
grid on
grid minor
box on
xlabel('X_2','fontsize',15); ylabel('R_X(X_2)')%,'rot',0)
set(gca,'fontsize',fs)
title('Nbhd 2')

subplot(2,2,3)
hold on
bar(xvec,bar2)
bar(xvec,bar1,'facecolor',.7*[1,1,1])
hold off
axis([.5,xvec(end)+.5,0,1])
set(gca,'XTick',[xvec(1),xvec(end)],'xticklabels',{'0','1'})
set(gca,'YTick',[0,1],'yticklabels',{'0','a'})
grid on
grid minor
box on
xlabel('X_1','fontsize',15); ylabel('R_X(X_1)')%,'rot',0)
set(gca,'fontsize',fs)

subplot(2,2,4)
hold on
bar(xvec,bar1)
bar(xvec,bar2,'facecolor',.7*[1,1,1])
hold off
axis([.5,xvec(end)+.5,0,1])
set(gca,'XTick',[xvec(1),xvec(end)],'xticklabels',{'0','1'})
set(gca,'YTick',[0,1],'yticklabels',{'0','a'})
grid on
grid minor
box on
xlabel('X_2','fontsize',15); ylabel('R_X(X_2)')%,'rot',0)
set(gca,'fontsize',fs)
%}