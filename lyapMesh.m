function f=lyapMesh
a=5:.25:7;
b=2:.5:4;
la=length(a); lb=length(b);
[x,y]=meshgrid(a,b);
out=zeros(la,lb,4);
for i=1:la
    ai=a(i);
    for j=1:lb
        bj=b(j);
        [f,g]=bnm2case2lyap(ai,bj);
        out(i,j,:)=g(end,:);
    end
end
fs=12; lw=2;
figure
for i=1:4
subplot(2,2,i)
surf(x,y,out(:,:,i)')
xlabel('\alpha')
ylabel('\beta')
zlabel(strcat('Lyapunov exponent',' ',num2str(i)))
axis tight
set(gca,'fontsize',fs);
grid on
grid minor
box on
end