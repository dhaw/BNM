function f=bnm2roots1b(a,g,b1,b2)
%Output is coefficients of p7 (symbolically)
%Symbolic substitution:
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);
B=xs*(xs-1);
[Q,R]=quorem(p,B,xs);
Q=simplify(Q);
R=simplify(R);
cq=coeffs(Q,xs);
Q=poly2sym(fliplr(cq),xs);%Fliplr - check ordering in commands/outputs
rq=roots(Q);
f=cq;
%
%a=1; g=2; b1=25; b2=50;
Qsub=subs(Q,[as,gs,b1s,b2s],[a,g,b1,b2]);
x=(0:.001:1);
z=subs(Qsub,x);
fs=20; lw=2; ms=10;
figure
hold on
plot(x,z,'k-','linewidth',1.5);
plot([0,1],[0,0],'k--','linewidth',1)
hold off
xlabel('X_1')
ylabel('p_7(X_1)')%,'rot',0)
set(gca,'fontsize',fs)
%Axis tight or axis([0,1,-1,1])
axis([0,1,-1,1])
grid on
grid minor
box on
%}
%{
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
%}
