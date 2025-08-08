function [P0, P1, P1min] = NullAndAlternativeDistributions_var_clust(b, bias, se, s, alpha, GeneralToSpecific)
%NULLANDALTERNATIVEDISTRIBUTIONS_VAR_CLUST
%   For cluster selection Monte Carlo experiments, returns:
%       P0    - null distribution of p-values (no p-hacking)
%       P1    - p-value distribution under thresholding (p-hacking)
%       P1min - p-value distribution under minimum approach (p-hacking)
%
%   Inputs:
%       b                 - [vector] regression coefficients
%       bias              - [matrix] vector of biases
%       se                - [matrix] vector of standard errors
%       s                 - 1 or 2, one- or two-sided test
%       alpha             - significance level
%       GeneralToSpecific - logical (1: general-to-specific, 0: specific-to-general)
%
%   Outputs:
%       P0, P1, P1min: [M x 1] vectors of p-value distributions for M specifications

    % Number of model specifications (columns of bias)
    M = size(bias, 2);

    % Initialize outputs
    P0    = zeros(M, 1);
    P1    = zeros(M, 1);
    
    % Compute t-statistics for all specifications
    t = (b + bias) ./ se;

    % Compute p-values depending on test type
    if s == 1
        % One-sided test
        p = 1 - normcdf(t);
    elseif s == 2
        % Two-sided test
        p = 2 * (1 - normcdf(abs(t)));
    else
        error('Input s must be 1 or 2');
    end

    % Null distribution (P0): no p-hacking
    if GeneralToSpecific == 1
        % Select the last model (most specific)
        P0 = p(end, :)';
    else
        % Select the first model (most general)
        P0 = p(1, :)';
    end

    % P1min: minimum p-value for each specification
    P1min = min(p)';

    % P1: p-value distribution under thresholding p-hacking
    for m = 1:M
        res = p(:, m);
        if GeneralToSpecific == 1
            res = flipud(res);
        end

        if min(res) > alpha
            P1(m) = min(res); % If no p-value passes alpha, take minimum
        else
            % Find the first p-value below alpha
            idx = find(res <= alpha, 1, 'first');
            P1(m) = res(idx);
        end
    end

end