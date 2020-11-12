function [fx,dfdx,dfdP] = f_onestatemodel(x,P,u,in)
As=(P(1));
Bs=(P(2));



fx=(As)*x(1)+Bs*(x(1)+u(1));
 dfdx =[];
 dfdP = [];