function f=Zplotevals(A)

eps=.01; thr=5;%Max beta+
al=(3+eps:eps:10)'; lal=length(al);
[bet1,bet2]=arrayfun(@betfromal,al);
c=find(al==4); bet1(c)=.5*(bet1(c-1)+bet1(c+1));%Lin interp b, not bet1?
thr=max(bet1)+thr;
bet2(bet2>thr)=thr; bet2(al>=4)=thr;
c=max(bet2); c=round(c/eps-1)*eps; bet=(9+eps:eps:c); lbet=length(bet);
%{
A=zeros(lal,lbet); B=A;
for i=1:lal%Are loops clonky here?
    ali=al(i);
    bet1i=bet1(i); bet1i=round(bet1i/eps)*eps; 
    beti=(bet1i:eps:bet2(i)); lbeti=length(beti);
    for j=1:lbeti-1
        betj=beti(j); b=round((betj-9)/eps);%+1
        try
            xij=roots([1,-2,(1+ali)/ali,(1-betj)/(ali*betj)]);
            xij=sortrows(xij); xij=xij(2); yij=xij*(1-xij);
            [d1]=arrayfun(@jac,xij,yij,ali,betj);
            A(i,b)=d1; %B(i,d)=d2;
        catch ME
            %A(i,b)=NaN;
        end
    end
end
A(A<=0)=NaN;
A(imag(A)~=0)=NaN;
f=A;%[al,bet1,bet2];
%}
fs=20; lw=2; ms=10;
[X,Y]=meshgrid(al,bet);
maxA=max(max(A));


x=(0:eps:10); lx=length(x); y=(0:eps:bet(end)); ly=length(y); AA=NaN(lx,ly); AA(end-lal+1:end,end-lbet+1:end)=A; A=AA;
figure
%surf(X,Y,A','EdgeColor','none')
%B=imagesc(al,bet,A');
hold on
B=imagesc(x,y,A');
set(B,'AlphaData',~isnan(A'))
plot([0,al(end)],[1,1],'k--','linewidth',lw)
plot(3,9,'ko','markersize',ms,'markerfacecolor','k','linewidth',1)
txt1='P';
text(2,9,txt1,'fontsize',fs)
xlabel('\alpha')
ylabel('\beta')
set(gca,'xtick',(0:2:10))
%set(gca,'ytick',[9,19,29,39])
set(gca,'FontSize',fs)
set(gca,'YDir','normal')
grid on
axis([0,al(end),0,bet(end)])%,0,maxA])
colorbar
caxis([0,maxA])
%}
end

function [f,g]=betfromal(al)
f=(9*al-2*al^2-2*sqrt(al*(al-3)^3))/(4-al);
g=(9*al-2*al^2+2*sqrt(al*(al-3)^3))/(4-al);
end

function [f]=jac(x,y,al,bet)
J=[x*(2-3*x)-y,-x;-y,bet*y*(2-3*al*y)-x];
eJ=eigs(J); eJ=sortrows(eJ);
%f=eJ(1); g=eJ(2);
f=det(J);
end