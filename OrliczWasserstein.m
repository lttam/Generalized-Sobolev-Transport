function val = OrliczWasserstein(phi, invphi, mu, nu, c, epsilon, tol)

if nargin < 7
    tol = 1e-6;
end

if nargin < 6
    epsilon = 0.1;
end

% phi: vectorize for a matrix
% invphi: only need for a scalar

dotp = @(x,y)sum(x(:).*y(:));
H = @(p)-sum( p(:).*(log(p(:)+1e-20)-1) );
options.niter = 1e+4; 
options.tau = 0;
options.verb = 0;

x_upp = max(c(:)) / invphi(1);

[~,~,gamma] = sinkhorn_log(mu,nu,c,epsilon,options);
ss = dotp(gamma, c);
x_low = (ss + (0.5*epsilon*(H(mu) + H(nu))))/invphi(1 + epsilon*(H(mu) + H(nu)));

cx_upp = phi(c ./ x_upp);
[~,~,gamma] = sinkhorn_log(mu,nu,cx_upp,epsilon,options);
fx_upp = dotp(gamma, cx_upp);

cx_low = phi(c ./ x_low);
[~,~,gamma] = sinkhorn_log(mu,nu,cx_low,epsilon,options);
fx_low = dotp(gamma, cx_low);

x_new = (x_low*fx_upp - x_upp*fx_low) / (fx_upp - fx_low);
iter = 1;
while (abs(x_upp - x_low) > tol)
    disp(['...' num2str(iter)]);
    if (x_new < x_upp) && (x_new > x_low)
        cx_new = phi(c ./ x_new);
        [~,~,gamma] = sinkhorn_log(mu,nu,cx_new,epsilon,options);
        fx_new = dotp(gamma, cx_new);

        if (fx_new < 1)
            x_upp = x_new;
            fx_upp = fx_new;
        else
            x_low = x_new;
            fx_low = fx_new;
        end
        x_new = (x_low*fx_upp - x_upp*fx_low) / (fx_upp - fx_low);
    else
        x_new = (x_low + x_upp)/2;
    end
    iter = iter + 1;
end

val = x_upp;

end



