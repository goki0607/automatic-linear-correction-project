function zero = newton_horner(f,x0,maxit,ftol,target)
% Horner's method augmented implementation of newton's method.
% Performs maxit iterations of newton's method on function
% f as defined below and uses function fdx to compute the 
% derivative of f at x. The function evaluations are done
% using an uncompensated Horner's algorithm.
% Inputs:      fx -- symbolic function to find the zero of
%              x0 -- initial starting point, defined by user
%           maxit -- the maximum number of iterations for the algorithm
%            ftol -- user specified tolerance
%          target -- known zero (we wish to find) for error analysis.
% Outputs:   zero -- the estimated solution of the zero of f.
f_x = expand(f);
f_dx = diff(f_x);
coef_x = sym2poly(f_x);
coef_dx = sym2poly(f_dx);
coef_x = single(coef_x);
coef_dx = single(coef_dx);
err = abs(x0-target);
sprintf('%.7e | %.7e | %.7e | %.7e | %.7e',0,x0,horner(coef_x,x0),...
    horner(coef_dx,x0),err)
for k=1:maxit
    f_x_res = horner(coef_x,x0);
    f_dx_res = horner(coef_dx, x0);
    x0 = x0 - f_x_res/f_dx_res;
    err = abs(x0-target);
    sprintf('%.7e | %.7e | %.7e | %.7e | %.7e',k,x0,f_x_res,f_dx_res,err)
    if f_x_res == 0
        disp('Found f(x) = 0!')
        break
    elseif abs(f_x_res) < ftol
        disp('|f| is within tolerance!')
        break
    end
end
zero=x0;
end 

function val = horner(vals,x0)
val = vals(1);
for i=2:length(vals)
    val = (val*x0)+vals(i);
end
end