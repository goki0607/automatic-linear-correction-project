function [A, b,x] = ill_matrix(n, alpha)
A = zeros(n,n);
A(1,1) = 100;
for i=2:n
    A(i,i) = 1;
end
for i=2:n
    for j=1:i-1
        A(i,j) = (-1)^(i+j)*alpha;
    end
end
b = zeros(n,1);
b(1) = 1;
for i=2:n
    b(i) = (-(alpha+1)/100)*(-2)^(i-2);
end
x = zeros(n,1);
x(1) = 0.01;
for i=2:n
    x(i) = (-1/100)*(-2)^(i-2);
end
        
            