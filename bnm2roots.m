function [f,g1]=bnm2roots
%Symbolic substitution:
syms xs zs as bs
p=subs(-as^2*(2*bs)*zs^3+as*(3*bs)*zs^2-(as+bs)*zs+1-xs,zs,-(2)*xs^3+(3)*xs^2-(1+1/as)*xs+1/as);
%syms xs zs as gs b1s b2s
%p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
c=coeffs(p,xs);
B=xs*(xs-1/2)*(xs-1);
[Q,R]=quorem(p,B,xs);
Q=simplify(Q);
R=simplify(R);
g1=subs(Q,xs,1/2);
cq=coeffs(Q,xs);
Q=poly2sym(fliplr(cq),xs);%Fliplr - check ordering in commands/outputs
rq=roots(Q);
f=cq;

cqdash=cq(2:end).*[1,2,3,4,5,6];%Derivative
p6dash=poly2sym(fliplr(cqdash),xs);
[p4,p4rem]=quorem(p6dash,xs-1/2,xs);
p4=simplify(p4);
p4rem=simplify(p4rem);
g1=coeffs(p4,xs);
%
a=9; g=1; b1=16; b2=b1;
Qsub=subs(Q,[as,bs],[a,b1]);
%Qsub=subs(Q,[as,gs,b1s,b2s],[a,g,b1,b2]);
x=(0:.001:1);
z=subs(Qsub,x);
fs=20; lw=2; ms=10;
figure
hold on
plot(x,z,'k-','linewidth',lw);
plot([0,1],[0,0],'k--','linewidth',1)
hold off
xlabel('X_1')
ylabel('p_6(X_1)')%,'rot',0)
set(gca,'fontsize',fs)
axis([0,1,-50,50])
grid on
grid minor
box on
%}
%{
syms xs zs as gs b1s b2s
p=subs(-as^2*(b1s+b2s)*zs^3+as*(b1s+2*b2s)*zs^2-(as+b2s)*zs+1-xs,zs,-(1+gs)*xs^3+(1+2*gs)*xs^2-(gs+1/as)*xs+1/as);
%}
