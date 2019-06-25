function [val, err] = two_sum(left, right)
% Error free transformation for addition/subtraction
% (OP) of two floating points numbers of any
% (but same) precision.
% Inputs:   left -- a floating point number such that
%                   left = fl(left)
%          right -- a floating point number such that
%                   right = fl(right).
% Let the error in fl(leftOPright) be expressed as:
% error = (leftOPright)-fl(leftOPright)
% Then the left of sum allows us to compute:
% (leftOPright) = fl(leftOPright) + error.
% Outputs:   val -- fl(leftOPright)
%            err -- error as defined above (i.e.
%                   elementary rounding error).
val = left+right;
s = val-left;
err = (left-(val-s))+(right-s);
end