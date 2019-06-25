function zero = newton_cena(f,x0,maxit,ftol,target,q)
% CENA implementation of newton's method 
% Performs maxit iterations of compensated newton's method
% on function f as defined below and compuets the derivative
% of x using Matlab's built in auto-differentiation package
% Inputs:       f -- symbolic function to find the zero of
%              x0 -- initial starting point, defined by user
%           maxit -- the maximum number of iterations for the algorithm
%            ftol -- user specified tolerance
%               q -- the # bits of the used precision's mantissa
%          target -- known zero (we wish to find) for error analysis.
% Outputs:   zero -- the estimated solution of the zero of f.
f_x = expand(f);
f_dx = diff(f_x);
coef_x = sym2poly(f_x);
coef_dx = sym2poly(f_dx);
coef_x = single(coef_x);
coef_dx = single(coef_dx);
for k=1:maxit
    [f_x_val, f_x_er] = comp_horner(coef_x, x0, q);
    [f_dx_val, f_dx_er] = comp_horner(coef_dx, x0, q);
    f_x_res = f_x_val+f_x_er;
    f_dx_res = f_dx_val+f_dx_er;
    [div, div_err] = approx_two_div(f_x_res, f_dx_res, q);
    div = div+div_err;
    [res, res_err] = two_sum(x0, -1*div);
    err = abs(x0-target);
    sprintf('%.7e | %.7e | %.7e | %.7e | %.7e',k,x0,f_x_res,f_dx_res,err)
    x0 = res+res_err;
    if f_x_res == 0
        disp('Found f(x) = 0!')
        break
    elseif abs(f_x_res) < ftol
        disp('|f| is within tolerance!')
        break
    end
end
zero = x0;
end

function [val, err] = comp_horner(vals, x, q)
val = vals(1);
pi = eye(1, length(vals)-1);
sigma = eye(1, length(vals)-1);
for i=2:length(vals)
    [mul_val, mul_err] = two_product(val, x, q);
    [sum_val, sum_err] = two_sum(vals(i), mul_val);
    val = sum_val;
    pi(i) = mul_err;
    sigma(i) = sum_err;
end
assert(length(pi)==length(sigma),'Issue at compensated Horner')
err = pi(length(pi))+sigma(length(sigma));
for i=2:length(pi)
    err = (pi(i)+sigma(i))+(err*x);
end
end