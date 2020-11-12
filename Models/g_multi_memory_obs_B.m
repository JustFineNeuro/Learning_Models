function [gx,dG_dX,dG_dPhi] = g_multi_memory_obs_B(Xt,Phi,ut,inG)


gx=-1*sign(ut(1,1)).*Xt(1); %%TODO: check sign flip is correct
gx(2)=Xt(2);
gx=sum(gx);


dG_dX = [];


if size(Phi,1) > 0
    dG_dPhi = zeros(size(Phi,1),size(G,1));
else
    dG_dPhi = [];
end