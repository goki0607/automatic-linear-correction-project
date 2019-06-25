function [val, err] = two_product(left, right, q)
% Error free transformation for multiplication (OP)
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
%                   elementary rounding error).
val = left*right;
[left_h, left_l] = split_high_low(left,q);
[right_h, right_l] = split_high_low(right,q);
err = left_l*right_l-(((val-left_h*right_h)-left_l*right_h)...
      -left_h*right_l);
end