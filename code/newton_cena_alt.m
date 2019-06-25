function zero = newton_cena_alt(f, x0, maxit, ftol, q)
f_x = expand(f);
f_dx = diff(f_x);
coef_x = sym2poly(f_x);
coef_dx = sym2poly(f_dx);
for k=1:maxit
    
    f_x_vals = eye(1, length(coef_x));
    f_x_err = 0;
   for i=1:length(coef_x)
       [pow_val, pow_err] = pow(x0, length(coef_x)-i, q);
       [mul_val, mul_err] = TwoProduct(coef_x(i), pow_val, q);
       f_x_vals(i) = mul_val;
       f_x_err = f_x_err + pow_err + mul_err;
   end
   [f_x_val, f_add_err] = sum(f_x_vals);
   f_x_err = f_x_err + f_add_err;
   
   f_dx_vals = eye(1, length(coef_dx));
   f_dx_err = 0;
   for i=1:length(coef_dx)
       [pow_val, pow_err] = pow(x0, length(coef_dx)-i, q);
       disp(pow_err)
       [mul_val, mul_err] = TwoProduct(coef_dx(i), pow_val, q);
       f_dx_vals(i) = mul_val;
       f_dx_err = f_dx_err + pow_err + mul_err;
   end
   [f_dx_val, f_add_err] = sum(f_dx_vals);
   f_dx_err = f_dx_err + f_add_err;
   
   f_x_true = f_x_val + f_x_err;
   f_dx_true = f_dx_val + f_dx_err;
   
   [div, div_err] = ApproxTwoDiv(f_x_true, f_dx_true, q);
   div = div+div_err;
   
   [t, e] = TwoSum(x0, -1*div);
   
   x0 = x0-f_x_true/f_dx_true
end
zero = x0;
end

function [val, err] = sum(vals)
val = 0;
err = 0;
for i=1:length(vals)
    [val, temp_err] = TwoSum(val, vals(i));
    err = err + temp_err;
end
end

function [val, err] = pow(num, n, q)
val = 1;
err = 0;
for i=1:n
    [val, temp_err] = TwoProduct(val, num, q);
    err = err + temp_err;
end
end