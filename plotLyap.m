function f=plotLyap(t,x)
fs=12; lw=2;
figure
hold on
plot(t,x(:,1),'-','linewidth',4)
plot(t,x(:,2),'--','linewidth',4)
plot(t,x(:,3),'-','linewidth',lw)
plot(t,x(:,4),'--','linewidth',lw)
xlabel('time')
ylabel('Lyapunov exponent')
axis tight
set(gca,'fontsize',fs);
grid on
grid minor
box on