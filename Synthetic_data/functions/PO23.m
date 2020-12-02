% l2/3_th - threshold function for L2/3 regularization
function [STq,pt,ss]=PO23(q,t)
pt=abs(q)-(2/3*t).^(3/4)*2;pt(pt>0)=1;pt(pt<=0)=0;
z=(1/16*q.^2+(q.^4/256-8*t.^3/729).^(1/2)).^(1/3)+(1/16*q.^2-(q.^4/256-8*t.^3/729).^(1/2)).^(1/3);
temp=((2*z).^(1/2)+(2*abs(q)./((2*z).^(1/2))-2*z).^(1/2)).^3;
STq=(1/8)*temp.*pt.*sign(q);
ss=abs(STq);