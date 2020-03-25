function [gx,dG_dX,dG_dPhi] = g_multi_memory_obs(Xt,Phi,ut,inG)


%%Here is where we would introdcued a gated representation between
%%competing system for parameter fitting%%%
%REquires gating input

idx=ut(2:3);
         gx=sum(Xt(idx));
%     end
   
% end


dG_dX = [];


if size(Phi,1) > 0
    dG_dPhi = zeros(size(Phi,1),size(G,1));
else
    dG_dPhi = [];
end