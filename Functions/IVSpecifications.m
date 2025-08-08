function [bias, se, Fstat] = IVSpecifications(eps, x, z, k)
%IVSPECIFICATIONS  All k-IV specifications: bias, SE, and first-stage F-stat
%
%   [bias, se, Fstat] = IVSpecifications(eps, x, z, k)
%
%   Inputs:
%       eps - [n x 1] regression residuals
%       x   - [n x 1] regressor of interest
%       z   - [n x K] available instruments
%       k   - number of instruments to include (k < K < 15)
%
%   Outputs:
%       bias  - [num_spec x 1] in-sample bias for each IV specification
%       se    - [num_spec x 1] in-sample standard error for each specification
%       Fstat - [num_spec x 1] first-stage F-statistic for each specification

    % Vector for selecting coefficient on x (second regressor)
    e2 = zeros(1, 2);
    e2(2) = 1;

    [n, K] = size(z);
    v = (1:K)';                               % Indices of available instruments
    specs = nchoosek(v, k);                   % All k-instrument combinations
    num_spec = size(specs, 1);

    bias   = zeros(num_spec, 1);
    se     = zeros(num_spec, 1);
    Fstat  = zeros(num_spec, 1);

    for i = 1:num_spec
        z_iv = [ones(n, 1), z(:, specs(i, :))];   % Selected IVs + intercept
        Pz = z_iv * ((z_iv' * z_iv) \ z_iv');     % Projection onto IV space

        X = [ones(n, 1), x];                     % Second stage: regressors (intercept + x)

        % IV bias for coefficient on x
        bias(i) = e2 * ((X' * Pz * X) \ (X' * Pz * eps));

        % IV regression residuals
        eps_hat_iv = X * ((X' * Pz * X) \ (X' * Pz * eps)) - eps;

        % IV standard error for coefficient on x
        se(i) = sqrt(e2 * ((inv(X' * Pz * X) * (sum((eps_hat_iv).^2) / (n - size(X, 2)))) * e2'));

        % First-stage: project x onto selected IVs
        pihat = (z_iv' * z_iv) \ (z_iv' * x);

        % Variance of first-stage coefficients
        V = (x' * (eye(n) - Pz) * x / (n - size(z_iv, 2))) * ((z_iv' * z_iv) \ eye(size(z_iv, 2)));

        % First-stage F-statistic (excludes intercept)
        Fstat(i) = pihat(2:end)' * (V(2:end, 2:end) \ pihat(2:end)) / (size(z_iv, 2) - 1);
    end

end
