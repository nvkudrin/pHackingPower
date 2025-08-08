function [k] = LagsBIC(y)
%LAGSBIC  Selects optimal Newey-West lag length via BIC
%
%   k = LagsBIC(y)
%
%   Inputs:
%       y - [T x 1] time series vector
%
%   Outputs:
%       k - scalar, optimal number of lags (integer between 0 and 4)
%
%   Description:
%       For lag lengths 0 to 4, computes the BIC for the corresponding Newey-West
%       autoregression and returns the lag length with minimum BIC.

    K = 4;                  % Maximum lag length considered
    T = length(y);          % Sample size
    bic = zeros(1, K+1);    % BIC for each lag (including zero)

    % Compute BIC for each lag (from 0 to 4)
    for j = 1:(K+1)
        Y = [];
        for i = 1:j
            % Construct lagged design matrix
            Y = [Y, y((j+1-i):(T-i+1))];
        end

        X = [ones(size(Y,1), 1), Y(:, 2:end)];  % Regressors: intercept + lags (excluding current y)
        M = eye(size(Y,1)) - X * ((X' * X) \ X');  % Projection matrix for residuals
        e = M * Y(:,1);                          % Regression residuals

        % BIC computation: log variance + penalty for number of parameters
        bic(j) = log(sum(e.^2) / length(e)) + (j+1) * log(T) / T;
    end

    % Return the lag length (0-based) minimizing BIC
    k = find(bic == min(bic)) - 1; % Adjust for MATLAB indexing (j = 1 is lag 0)
end
