function x = substitution_cena(A, b, q)
[m,n] = size(A);
x = zeros(m, 1);
x_corr = zeros(m, 1);
for i=1:n
    sum = b(i);
    u = 0;
    for j=1:i-1
        [mul, mul_err] = two_product(A(i,j), x(j), 2);
        [sum, sum_err] = two_sum(sum, -mul);
        u = u-(A(i,i)*x_corr(j)+mul_err-sum_err);
    end
    [x(i), x_err_i] = approx_two_div(sum, A(i,i), q);
    x_corr(i) = u/A(i,i)+x_err_i;
end
disp(x)
disp(x_corr)
for i=2:n
   x(i) = x(i)+x_corr(i);
end
end
