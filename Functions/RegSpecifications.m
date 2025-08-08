function [bias, se] = RegSpecifications(eps, x, z, k)
%REGSPECIFICATIONS Computes in-sample bias and standard error for all combinations of controls.
%   [bias, se] = RegSpecifications(eps, x, z, k)
%
%   Inputs:
%       eps : [n x 1] vector of regression residuals
%       x   : [n x 1] regressor of interest
%       z   : [n x K] matrix of additional controls
%       k   : scalar, number of controls to include (k < K < 15)
%
%   Outputs:
%       bias : [num_spec x 1] vector of in-sample bias for each specification
%       se   : [num_spec x 1] vector of standard errors for each specification
%
%   Examines all combinations of k controls and returns bias and se for the coefficient on x.

    % Set up selector vector for coefficient on x
    e1 = zeros(1, k + 2);
    e1(1) = 1; % Picks out the coefficient on x

    [n, K] = size(z);

    if k > 0
        % Get all combinations of k controls out of K available
        v = (1:K)';
        specs = nchoosek(v, k);       % [num_spec x k] indices of selected controls
        num_spec = size(specs, 1);

        bias = zeros(num_spec, 1);
        se   = zeros(num_spec, 1);

        for i = 1:num_spec
            % Construct design matrix: [x, intercept, selected controls]
            X = [x, ones(n, 1), z(:, specs(i, :))];
            % Compute bias for the coefficient on x
            bias(i) = e1 * ((X' * X) \ (X' * eps));
            % Calculate regression residuals
            eps_hat = X * ((X' * X) \ (X' * eps)) - eps;
            % Compute standard error for coefficient on x
            se(i) = sqrt(e1 * (inv(X' * X) * (sum(eps_hat.^2) / (n - size(X, 2)))) * e1');
        end

    else
        % No controls: only x and intercept
        X = [x, ones(n, 1)];
        bias = e1 * ((X' * X) \ (X' * eps));
        eps_hat = X * ((X' * X) \ (X' * eps)) - eps;
        se = sqrt(e1 * (inv(X' * X) * (sum(eps_hat.^2) / (n - size(X, 2)))) * e1');
    end

end