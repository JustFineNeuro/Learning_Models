function [fx,dfdx,dfdP] = f_twomemorymodel_B(x,P,u,in)
As=sig_trans(P(1));
Af=sig_trans(P(2));
Bs=-sig_trans(P(3));
Bf=-sig_trans(P(4));
% As=(P(1));
% Af=(P(2));
% Bs=(P(3));
% Bf=(P(4));
%REquires gating input


pe=u(1)-sum(x);



fx=[(As) 0; 0 (Af) ]*[x(1);x(2)]+[(Bs)*pe;(Bf)*pe];
dfdx = [];
dfdP = [];