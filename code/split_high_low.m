function [p1, p2] = split_high_low(val, q)
% Split function that splits a floating point number
% to two parts, such that p1+p2=val with each part
% having at most r-1 non-zero bits, r is defined as
% the floor of q (i.e. the precision of val).
% Inputs:   val -- sum number val = fl(val)
%             q -- the # bits of val's mantissa.
% Outputs:   p1 -- the first part of val
%            p2 -- the second part of val.
r = ceil(q/2);
inter = val*((2^r)+1);
p1 = inter-(inter-val);
p2 = val-p1;
end