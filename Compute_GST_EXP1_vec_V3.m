function val_vec = Compute_GST_EXP1_vec_V3(XX, WW)

% XX: k x dim
% WW: k x dim
% val_vec: k x 1

constEPS = 1e-10;

kk = size(XX, 1);
val_vec = zeros(kk, 1);

options = optimoptions('fmincon', 'Algorithm','trust-region-reflective',...
    'FunctionTolerance', 1e-3, ...
    'ConstraintTolerance', 1e-3, ...
    'MaxFunctionEvaluations', 100, ...
    'MaxIterations', 100, ...
    'OptimalityTolerance', 1e-3, ...
    'StepTolerance', 1e-3, ...
    'UseParallel', true, ...
    'SpecifyObjectiveGradient',true, ...
    'HessianFcn','objective', ...
    'Display', 'off');
k0 = 1;

for ii = 1:kk
    xii = XX(ii, :);
    wii = WW(ii, :);

    if sum(abs(xii)) < constEPS   
        val_vec(ii) = 0;
    else
        % solve 1d opt
        func = @(k) obj_func_Exp1(k, wii', xii');                
        k_start = fmincon(func, k0,[],[],[],[],constEPS,inf, [], options);

        val_vec(ii) = obj_func_Exp1(k_start, wii', xii');
    end
end

end

% compute objective function & gradient function
function [f, g, h] = obj_func_Exp1(k, w, x)
    % obj
    kx = k*x;
    phi_kx = exp(kx) - kx - 1;
    wphi = w' * phi_kx;
    f = (1.0 + wphi)/k;
    % grad
    if nargout > 1
        g1 = -(1.0 + wphi)/(k^2);

        gradphi_kx = exp(kx) - 1;
        wxgradphi_kx = w' * (x .* gradphi_kx);
        g = g1 + (wxgradphi_kx/k);

        if nargout > 2

            h1 = 2*(1+wphi)/(k^3);
            h2 = -2*wxgradphi_kx/(k^2);

            hessphi_kx = exp(kx);
            h = h1 + h2 + ((w' * (x .* x .* hessphi_kx))/k);

        end
    end
end

