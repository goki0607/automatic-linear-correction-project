function [val, err] = approx_two_div(left, right, q)
% NON-error free transformation for division (OP)
% (since an error free transformation is not possible)
% of two floating points numbers of any (but same)
% precision.
% Inputs:   left -- a floating points number such that
%                   left = fl(left)
%          right -- a floating point number such that
%                   right = fl(right)
%              q -- the # bits of left and rights's
%                   mantissas this implies both left 
%                   and right must be of the same 
%                   precision.
% Outputs:   val -- fl(leftOPright)
%            err -- error as defined above (i.e.
%                   approx. division rounding error).
assert(right~=0,'A division by zero has occurred!')
val = left/right;
[temp_val, temp_err] = two_product(val,right,q);
err = (temp_val-left-temp_err)/right;
end
