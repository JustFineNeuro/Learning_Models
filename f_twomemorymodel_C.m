function [fx,dfdx,dfdP] = f_twomemorymodel_C(x,P,u,in)
As=sig_trans(P(1));
Af=sig_trans(P(2));
Bs=-sig_trans(P(3));
Bf=-sig_trans(P(4));
Gen=sig_trans(P(5));
% As=(P(1));
% Af=(P(2));
% Bs=(P(3));
% Bf=(P(4));
%REquires gating input
% persistent t
% if in.active(:,t)==1
%     idx=[1 3];
% else
%     idx=[2 3];
% end
idx=u(2:3);
pe=u(1)-sum(x(idx));
Ax=[(As) 0 0; 0 As 0; 0 0 (Af) ]*[x(1);x(2);x(3)];

BF_pe=(Bf)*pe;

if idx(1)==1

    BS_pe=(Bs)*pe;
    fx=Ax+[BS_pe;0;BF_pe];
else
    BS_pe=(Bs)*pe;
    fx=Ax+[0;BS_pe;BF_pe];
end

% fx=[(As) 0 0; 0 As 0; 0 0 (Af) ]*[x(1);x(2);x(3)]+[(Bs)*pe;(Bs)*pe;(Bf)*pe];
dfdx = [];
dfdP = [];