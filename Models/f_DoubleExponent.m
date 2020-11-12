function [fx] = f_DoubleExponent(x,P,u,in)
SpanSlow=(P(1));
SpanFast=(P(2));
KSlow=(P(3));
KFast=(P(4));
Plateau=P(5);

fx=Plateau + SpanFast*exp(-KFast*x(1)) + SpanSlow*exp(-KSlow*x(1));

end

