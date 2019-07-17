function [f,g,ic1,ic2]=bnm2case2lyapABfixed(a,b)

runs=20;
lyap=zeros(runs,4);
out=zeros(runs,6);%x1z1x2z2x3z3
xic=zeros(runs,3);
zic=xic;
for i=1:runs
    x0=rand(1,3); x0=x0/sum(x0);
    z0=rand(1,3); z0=z0/sum(z0)/a;
    xic(i,:)=x0;
    zic(i,:)=z0;
    [~,g]=bnm2case2lyap(a,b,x0,z0);%Set to rand ICs
    lyap(i,:)=g(end,:);
    f=bnm2case2diffTS(a,b,x0,z0,0);
    out(i,:)=f(end,:);
end
f=lyap;
g=out;
ic1=xic;
ic2=zic;

lyap2=lyap;
%lyap2(lyap2<-.1)=-.1;
%lyap2(lyap2>.1)=.1;
runvec=1:runs;
cmap1=zeros(4,3);
cmap2=[0,0,1;1,0,0;0,0,1;1,0,0;0,0,1;1,0,0];

figure
fs=12; lw=2;
subplot(2,1,1)
hold on
plot([0,runs+1],[0,0],'k-','linewidth',1.5);
h1=scatter(runvec,lyap2(:,1),[],cmap1(1,:),'o','linewidth',lw);%'filled','o'
h2=scatter(runvec,lyap2(:,2),[],cmap1(2,:),'x','linewidth',lw);
h3=scatter(runvec,lyap2(:,3),[],cmap1(3,:),'^','linewidth',lw);
h4=scatter(runvec,lyap2(:,4),[],cmap1(4,:),'s','linewidth',lw);
hold off
xlabel('Run');
ylabel('Lyap. exp.');
set(gca,'FontSize',fs);
axis tight%([0,runs+1,-.1,.1])
legend([h1,h2,h3,h4],'\lambda_1','\lambda_2','\lambda3','\lambda_4','location','northeastoutside')
grid on
grid minor
box on

subplot(2,1,2)
hold on
h1=scatter(runvec,out(:,1),[],cmap2(1,:),'filled','o','linewidth',lw);%'filled','o'
h2=scatter(runvec,a*out(:,2),[],cmap2(2,:),'filled','o','linewidth',lw);
h3=scatter(runvec,out(:,3),[],cmap2(3,:),'o','linewidth',lw);
h4=scatter(runvec,a*out(:,4),[],cmap2(4,:),'o','linewidth',lw);
h5=scatter(runvec,out(:,5),[],cmap2(5,:),'x','linewidth',lw);
h6=scatter(runvec,a*out(:,6),[],cmap2(6,:),'x','linewidth',lw);
hold off
xlabel('Run');
ylabel('Population');
set(gca,'FontSize',fs);
axis ([0,runs+1,0,1])
legend([h1,h2,h3,h4,h5,h6],'X_1','Z_1','X_2','Z_2','X_3','Z_3','location','northeastoutside')
grid on
grid minor
box on