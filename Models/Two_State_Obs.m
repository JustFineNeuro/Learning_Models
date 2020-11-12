function [gx,dG_dX,dG_dPhi] = Two_State_Obs(Xt,Phi,ut,inG)

% if exist('inG')
%     if ~isempty(inG.statemap)
%         G=transpose(inG.statemap(:));
%     else
         G=ones(1,size(Xt,1));
% %     end
%     
% %     try
% %   'worked'
% %         gx = (G*Xt)+ut(1);
% %     catch
%         
         gx=G*Xt;
%     end
   
% end


dG_dX = [];


if size(Phi,1) > 0
    dG_dPhi = zeros(size(Phi,1),size(G,1));
else
    dG_dPhi = [];
end