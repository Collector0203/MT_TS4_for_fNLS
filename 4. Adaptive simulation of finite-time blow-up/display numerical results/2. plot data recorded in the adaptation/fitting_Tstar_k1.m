function result = fitting_Tstar_k1(T_fix, U_fix)
%FITTING_TSTAR_K1
% Fit the model:
%   log(U_fix) = -k1 * log(Tstar - T_fix) + k2
%
% Inputs:
%   T_fix : Nx1 or 1xN array
%   U_fix : Nx1 or 1xN array, must be > 0
%
% Output:
%   result is a struct with fields:
%       result.Tstar
%       result.k1
%       result.k2
%       result.t
%       result.u
%       result.y
%       result.yfit
%       result.Ufit
%       result.res_log
%       result.SSE_log
%       result.RMSE_log
%       result.R2_log
%       result.method
%       result.exitflag
%       result.output
%       result.valid_idx
%       result.bestResidual
%
% Strategy:
%   - Optimize only Tstar and k1 numerically
%   - Compute k2 analytically for each (Tstar, k1)
%   - Enforce Tstar > max(T_fix) by parameterization:
%         Tstar = max(T_fix) + exp(alpha)
%   - Enforce k1 > 0 by parameterization:
%         k1 = exp(beta)

    %% Input preprocessing
    if nargin < 2
        error('Both T_fix and U_fix must be provided.');
    end

    t = double(T_fix(:));
    u = double(U_fix(:));

    valid_idx = isfinite(t) & isfinite(u) & (u > 0);
    t = t(valid_idx);
    u = u(valid_idx);

    if isempty(t)
        error('No valid data found. U_fix must be positive and finite.');
    end

    [t, order] = sort(t);
    u = u(order);

    y = log(u);
    tmax = max(t);
    trng = max(t) - min(t);
    scale = max([1; abs(t)]);

    %% Parameterization setup
    % Tstar = tmax + exp(alpha)
    dmin = max(1e-14 * scale, 10 * eps(tmax));
    dmax = max(10 * max(trng, 1), 1e6 * dmin);

    % k1 = exp(beta). The interval below is only used to generate
    % coarse initial guesses; the subsequent local optimization is unbounded
    % in beta.
    k1_min = 1e-4;
    k1_max = 1e4;

    %% Coarse scan for robust initial guesses
    alpha_list = linspace(log(dmin), log(dmax), 220);
    beta_list  = linspace(log(k1_min), log(k1_max), 220);

    SSE_grid = inf(numel(alpha_list), numel(beta_list));

    for i = 1:numel(alpha_list)
        alpha = alpha_list(i);
        Tstar_tmp = tmax + exp(alpha);

        z = log(Tstar_tmp - t);

        for j = 1:numel(beta_list)
            beta = beta_list(j);
            k1_tmp = exp(beta);

            k2_tmp = mean(y + k1_tmp * z);
            y_tmp = -k1_tmp * z + k2_tmp;
            r_tmp = y - y_tmp;

            SSE_grid(i, j) = sum(r_tmp.^2);
        end
    end

    valid_mask = isfinite(SSE_grid);
    if ~any(valid_mask(:))
        error('No valid coarse-scan points were found.');
    end

    [~, sorted_idx] = sort(SSE_grid(:), 'ascend');
    nStarts = min(12, nnz(valid_mask));
    start_idx = sorted_idx(1:nStarts);

    X0 = zeros(nStarts, 2);
    for k = 1:nStarts
        [ia, ib] = ind2sub(size(SSE_grid), start_idx(k));
        X0(k, 1) = alpha_list(ia);
        X0(k, 2) = beta_list(ib);
    end

    %% Multi-start local refinement
    use_lsqnonlin = (exist('lsqnonlin', 'file') == 2);

    bestSSE = inf;
    bestx = [NaN; NaN];
    bestResidual = [];
    bestExitflag = NaN;
    bestOutput = [];
    bestMethod = '';

    for k = 1:nStarts
        x0 = X0(k, :).';

        if use_lsqnonlin
            resfun = @(x) residual_fun(x, t, y, tmax);

            opts = optimoptions('lsqnonlin', ...
                'Display', 'off', ...
                'FunctionTolerance', 1e-14, ...
                'StepTolerance', 1e-14, ...
                'OptimalityTolerance', 1e-14, ...
                'MaxIterations', 2e4, ...
                'MaxFunctionEvaluations', 2e5);

            [xopt, ~, residual, exitflag, output] = lsqnonlin(resfun, x0, [], [], opts);
            SSE = sum(residual.^2);
            method = 'lsqnonlin';
        else
            objfun = @(x) scalar_obj(x, t, y, tmax);

            opts = optimset( ...
                'Display', 'off', ...
                'TolX', 1e-14, ...
                'TolFun', 1e-14, ...
                'MaxIter', 2e4, ...
                'MaxFunEvals', 2e5);

            [xopt, SSE, exitflag, output] = fminsearch(objfun, x0, opts);
            residual = residual_fun(xopt, t, y, tmax);
            method = 'fminsearch';
        end

        if isfinite(SSE) && (SSE < bestSSE)
            bestSSE = SSE;
            bestx = xopt;
            bestResidual = residual;
            bestExitflag = exitflag;
            bestOutput = output;
            bestMethod = method;
        end
    end

    if ~all(isfinite(bestx))
        error('Optimization failed to find a valid solution.');
    end

    %% Recover fitted parameters
    alpha = bestx(1);
    beta = bestx(2);

    Tstar = tmax + exp(alpha);
    k1_fit = exp(beta);

    z = log(Tstar - t);
    k2_fit = mean(y + k1_fit * z);

    yfit = -k1_fit * z + k2_fit;
    Ufit = exp(yfit);
    res_log = y - yfit;

    SSE = sum(res_log.^2);
    RMSE_log = sqrt(mean(res_log.^2));

    denom = sum((y - mean(y)).^2);
    if denom > 0
        R2_log = 1 - SSE / denom;
    else
        R2_log = NaN;
    end

    %% Build output struct
    result = struct();
    result.Tstar = Tstar;
    result.k1 = k1_fit;
    result.k2 = k2_fit;

    result.t = t;
    result.u = u;
    result.y = y;
    result.yfit = yfit;
    result.Ufit = Ufit;
    result.res_log = res_log;

    result.SSE_log = SSE;
    result.RMSE_log = RMSE_log;
    result.R2_log = R2_log;

    result.method = bestMethod;
    result.exitflag = bestExitflag;
    result.output = bestOutput;

    result.valid_idx = valid_idx;
    result.bestResidual = bestResidual;
end

function r = residual_fun(x, t, y, tmax)
% Return the residual vector in log-domain

    alpha = x(1);
    beta = x(2);

    Tstar = tmax + exp(alpha);
    k1 = exp(beta);

    if ~isfinite(Tstar) || ~isfinite(k1)
        r = 1e6 * ones(size(y));
        return;
    end

    d = Tstar - t;
    if any(d <= 0) || any(~isfinite(d))
        r = 1e6 * ones(size(y));
        return;
    end

    z = log(d);
    if any(~isfinite(z))
        r = 1e6 * ones(size(y));
        return;
    end

    k2 = mean(y + k1 * z);
    ymodel = -k1 * z + k2;

    r = y - ymodel;
    if any(~isfinite(r))
        r = 1e6 * ones(size(y));
    end
end

function J = scalar_obj(x, t, y, tmax)
% Scalar objective for fminsearch

    r = residual_fun(x, t, y, tmax);
    J = sum(r.^2);
end