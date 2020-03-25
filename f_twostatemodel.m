function [fx,dfdx,dfdP] = f_twostatemodel(x,P,u,in)
As=(P(1));
Af=P(2);
Bs=(P(3));
Bf=P(4);

%REquires gating input
idx=in.active;
pe=sum(x(idx))-u(1);

fx=[(As) 0 0; 0 As 0; 0 0 (Af) ]*[x(1);x(2);x(3)]+[(Bs)*pe;(Bs)*pe;(Bf)*pe];
dfdx = [];
dfdP = [];