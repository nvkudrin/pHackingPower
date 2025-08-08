function [P0, P1, P1min] = NullAndAlternativeDistributions_var_laglength(b, bias, se, s, alpha, bic)
%NULLANDALTERNATIVEDISTRIBUTIONS_VAR_LAGLENGTH
%   For lag length selection Monte Carlo experiments, returns:
%       P0    - null distribution of p-values (no p-hacking, BIC selected lag)
%       P1    - p-value distribution under thresholding (p-hacking)
%       P1min - p-value distribution under minimum approach (p-hacking)
%
%   Inputs:
%       b     - [vector] regression coefficients
%       bias  - [matrix] bias values
%       se    - [matrix] standard errors
%       s     - 1 or 2 (one- or two-sided test)
%       alpha - significance level
%       bic   - [vector] BIC-selected lag indices (integer, 0-based, length M)
%
%   Outputs:
%       P0, P1, P1min: [M x 1] vectors of p-value distributions for M specifications

    M = size(bias, 2);        % Number of specifications (columns of bias)
    P0 = zeros(M, 1);         % Null p-value distribution
    P1 = zeros(M, 1);         % Thresholding p-hacking p-value distribution

    % Compute t-statistics for all specifications
    t = (b + bias) ./ se;

    % Compute p-values according to test type
    if s == 1
        % One-sided test
        p = 1 - normcdf(t);
    elseif s == 2
        % Two-sided test
        p = 2 * (1 - normcdf(abs(t)));
    else
        error('Input s must be 1 or 2.');
    end

    bic = bic + 1; % to create valid Matlab indices

    % Null distribution (P0): use p-value from BIC-selected model for each m
    for m = 1:M
        P0(m) = p(bic(m), m); % Note: bic(m) should be a valid row index for p(:, m)
    end

    % P1min: minimum p-value for each specification
    P1min = min(p)';

    % P1: thresholding p-hacking (search over all models from BIC-selected lag onwards)
    for m = 1:M
        res = p(bic(m):end, m);

        if min(res) > alpha
            P1(m) = min(res); % If no p-value passes alpha, take the minimum
        else
            % Find the first p-value below alpha
            idx = find(res <= alpha, 1, 'first');
            P1(m) = res(idx);
        end
    end

end
