function [x_sig] = sig_trans(x_native)
x_sig=1./(1+exp(-x_native));
end

