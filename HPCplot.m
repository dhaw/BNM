function f=HPCplot(f1,f2)
fs=12; lw=2;
figure
%For UN (f1,f2)

M=f1;
a=(2:2:30);
n=(2:2:10);
la=length(a);
ln=length(n);
lm=size(M,3);%Input only #layers need
cmap=.5*colormap(jet(ln));
hold all
for k=1:lm
    for j=2:ln
        plot(a,M(:,j,k),'ko');%,'MarkerFaceColor','k')%'markersize',2)%,'LineWidth',1
    end
end
h1=plot(a,f2(:,1),'o-','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
for j=2:ln-1
    h3=plot(a,f2(:,j),'o-','color',cmap(j,:),'LineWidth',3,'MarkerFaceColor',cmap(j,:));
end
h2=plot(a,f2(:,end),'o-','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
hold off
axis([0,max(a),0,1.05])
xlabel('a','FontSize',15)
ylabel('Dissimilarity','FontSize',fs)
set(gca,'FontSize',fs)
legend([h1 h2],{'n=2','n=10'},'location','W');
grid on
grid minor
box on
%}

%For RL (f,g,h,xs,ys)
%{
n=(4:2:30);
m=(1:5);
ln=length(n);
lm=length(m);
cmap=.5*colormap(jet(lm));
%{
hold on
for i=1:lm
    plot(n,f(:,i),'x:','color',cmap(i,:));
end
plot([1,n(end)],[0,0],'k:','LineWidth',2)
xlabel('n','FontSize',15)
ylabel('Dissimilarity','FontSize',15)
set(gca,'FontSize',15)
axis([0,n(end),0,1])
hold off
%}
hold on
plot([1,n(end)],[0,0],'k:','LineWidth',2)
h1=plot(n,g(:,1),'o:','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
h3=plot(n,h(:,1),'o-','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
for i=2:lm-1
    plot(n,g(:,i),'o:','color',cmap(i,:),'LineWidth',3,'MarkerFaceColor',cmap(i,:));
    plot(n,h(:,i),'o-','color',cmap(i,:),'LineWidth',3,'MarkerFaceColor',cmap(i,:));
end
h2=plot(n,g(:,end),'o:','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
h4=plot(n,h(:,end),'o-','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
xlabel('n','FontSize',15)
ylabel('M,M^{XY}','FontSize',15)
set(gca,'FontSize',15)
axis([0,n(end),min(g(:,1)),max(h(:,1))])
legend([h3 h4],{'m=1','m=5'},'location','NW');
hold off
%}

%For SF/star (f,g,h,xs,ys,in1,xyin,xh1)
%{
f=nanmean(f,2);
g=nanmean(g,2);
h=nanmean(h,2);
%xin=nanmean(xin,2);
%xyin=nanmean(xyin,2);
in1=nanmean(in1,2);
a=(2:1:30);
%a=(2:2:20);
in1=in1./a';
la=length(a);
%cmap=.5*colormap(jet(2));
%hold on
plot(a,f,'ko-','LineWidth',3,'MarkerFaceColor','k')
xlabel('a','FontSize',15)
ylabel('Dissimilarity','FontSize',15)
set(gca,'FontSize',15)
axis ([0,a(end),0,1])
%hold off
figure
hold on
plot([0,a(end)],[0,0],'k:','LineWidth',2);
plot(a,g,'ko-','LineWidth',3,'MarkerFaceColor','k')
plot(a,h,'ko:','LineWidth',3,'MarkerFaceColor','k')
plot(a,nanmean(in1,2),'o:','color',[0,.5,0],'LineWidth',3,'MarkerFaceColor',[0,.5,0]) %in1
xlabel('a','FontSize',15)
ylabel('M,M^{XY}','FontSize',15)
set(gca,'FontSize',15)
mm=max(max(h),max(in1));
axis tight%([0,a(end),min(g)-.5,mm])
legend('M','M^{XY}','P_1','location','SW')
hold off
%{
figure
hold on
plot(a,xyin,'ko-','LineWidth',2,'MarkerFaceColor','k')
plot(a,xin,'ko:','LineWidth',2,'MarkerFaceColor','k')
xlabel('n','FontSize',15)
ylabel('Number of Occupied Cells','FontSize',15)
set(gca,'FontSize',15)
axis([0,a(end),0,max(xyin)])
legend('X&Y','X only')
hold off
%}
%}

%For RG/SW (f,g,h,xs,ys)
%{
n=50;
a=(2:2:20);%********
p=(0:.05:1);%(.1:.1:1)(0:.05:1);%******** %start at 0 for SW, .1 for RG
la=length(a);
lp=length(p);
cmap=.5*colormap(jet(la));

hold on
h1=plot(p,f(1,:),'o-','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
for i=2:la-1
    plot(p,f(i,:),'o-','color',cmap(i,:),'LineWidth',3,'MarkerFaceColor',cmap(i,:))
end
h2=plot(p,f(end,:),'o-','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
plot([1,n(end)],[0,0],'k:','LineWidth',2)
xlabel('p','FontSize',15)
ylabel('Dissimilarity','FontSize',15)
set(gca,'FontSize',15)
axis([0,1,0,1.01])
legend([h1 h2],{'a=2','a=20'},'location','SW');
hold off

figure
hold on
plot([0,p(end)],[0,0],'k:','LineWidth',2)
h1=plot(p,g(1,:),'o:','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
h3=plot(p,h(1,:),'o-','color',cmap(1,:),'LineWidth',3,'MarkerFaceColor',cmap(1,:));
for i=1:la
    plot(p,g(i,:),'o:','color',cmap(i,:),'LineWidth',3,'MarkerFaceColor',cmap(i,:))
    plot(p,h(i,:),'o-','color',cmap(i,:),'LineWidth',3,'MarkerFaceColor',cmap(i,:))
end
h2=plot(p,g(end,:),'o:','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
h4=plot(p,h(end,:),'o-','color',cmap(end,:),'LineWidth',3,'MarkerFaceColor',cmap(end,:));
xlabel('p','FontSize',15)
ylabel('M,M^{XY}','FontSize',15)
set(gca,'FontSize',15)
ming=min(min(g));
maxh=max(max(h));
axis([0,1,ming,maxh+.1])
legend([h3 h4],{'a=2','a=20'},'location','E');
hold off
%}

%For control (f)
%{
%f=nanmean(f,3);
n=(2:2:20);%********
c=(1:5);%********
ln=length(n);
lc=length(c);
cmap=.5*colormap(jet(lc));
hold on
for i=1:lc
    plot(n,f(:,i),'o-','color',cmap(i,:),'LineWidth',3)
end
%plot([1,n(end)],[0,0],'k:','LineWidth',2)
xlabel('n','FontSize',15)
ylabel('Dissimilarity','FontSize',15)
set(gca,'FontSize',15)
axis([0,max(n),0,1.1])
legend('c=1','c=n','c=2n','c=3n','c=n^2','location','SE')%'northeastoutside')
hold off
%}