function zero = newton_naive(fx,fdx,x0,maxit,ftol,target)
% Naive implementation of newton's method.
% Performs maxit iterations of newton's method on function
% f as defined below and uses function fdx to compute the 
% derivative of f at x.
% Inputs:      fx -- anonymous function to find the zero of
%             fdx -- derivative of f (also anonymous)
%              x0 -- initial starting point, defined by user
%           maxit -- the maximum number of iterations for the algorithm
%            ftol -- user specified tolerance
%          target -- known zero (we wish to find) for error analysis.
% Outputs:   zero -- the estimated solution of the zero of f.
err = abs(x0-target);
sprintf('%.7e | %.7e | %.7e | %.7e | %.7e',0,x0,fx(x0),fdx(x0),err)
for k=1:maxit
    x0 = x0 - fx(x0)/fdx(x0);
    err = abs(x0-target);
    sprintf('%.7e | %.7e | %.7e | %.7e | %.7e',k,x0,fx(x0),fdx(x0),err)
    if fx(x0) == 0
        disp('Found f(x) = 0!')
        break
    elseif abs(fx(x0)) < ftol
        disp('|f| is within tolerance!')
        break
    end
end
zero=x0;
end   