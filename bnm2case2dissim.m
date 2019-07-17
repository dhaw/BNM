function f=bnm2case2dissim%(x0,y0)
x0=rand(1,3); x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
y0=rand(1,3); y0=y0/sum(y0); %[z1,z2,z3] - can input as ratio

thismany=20;
%If 2x1nbhd:
%{
x0=x0([1,2]); y0=y0([1,2]);
%}
%
%IC mesh:
xv=0:.05:1; lx=length(x);
yx=x;
[xm,ym]=meshgrid(xv,yv);
xx=reshape(xm,lx^2,1);
yy=reshape(ym,lx^2,1);
ic=[xx,yy];
%}
a=.1:.1:10;
la=length(a);
b=.2:.2:40;
lb=length(b);
[A,B]=meshgrid(a,b);
%X2=1-X; Y2=1/a-Y;
%Xsum=X+X2; Ysum=Y+Y2
D=zeros(la,lb,thismany);
for i=1:la
    ai=a(i);
    y0a=y0/ai;
    %
    icab=ic;
    icab(icab(:,2)>1/ai)=[];
    thismany=length(icab);
    %}
    for j=1:lb
        Dijk=zeros(thismany,1);
        parfor k=1:thismany
            %{
            x0=rand(1,3); x0=x0/sum(x0); %[x1,x2,x3], normalised just in case
            y0=rand(1,3); y0a=y0/sum(y0)/ai; %[z1,z2,z3] - can input as ratio
            %If 2x1nbhd:
            %
            x0=x0([1,2]); y0a=y0a([1,2]);
            %}
            %}
            dij=bnm2case2diffTSevents(ai,b(j),x0,y0a,0);%x1y1x2y2
            %dij=bnm2case2diffTSevents2x1d(ai,b(j),x0,y0a,0);%x1y1x2y2
            Dijk(k)=(abs(dij(1)-ai*dij(2))+abs(dij(3)-ai*dij(4)))/2;
        end
        D(i,j,:)=Dijk;
    end
end
D=nanmean(D,3);
f=D;%[x0;y0a];
%
fs=20; lw=2; ms=10;
figure
%h1=imagesc(a,b,D');
h1=contourf(a,b,D','LineStyle','none');
%set(D,'AlphaData',~isnan(D'))
%plot([0,a(end)],[1,1],'k--','linewidth',lw)
%plot(3,9,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
%txt1='P';
%text(2,9,txt1,'fontsize',fs)
xlabel('\alpha')
ylabel('\beta','rot',0)
%set(gca,'xtick',(0:2:10))
%set(gca,'ytick',[9,19,29,39])
set(gca,'FontSize',fs)
set(gca,'YDir','normal')
grid on
axis([0,a(end),0,b(end)])%,0,maxA])
colorbar
caxis([0,1])
