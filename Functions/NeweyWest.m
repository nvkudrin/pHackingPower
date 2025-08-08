function [Q] = NeweyWest(eps, x, L)
%NEWEYWEST Calculates Newey-West standard error for a regressor x
%
%   Q = NeweyWest(eps, x, L)
%
%   Inputs:
%       eps - [T x 1] vector of regression residuals
%       x   - [T x 1] regressor
%       L   - integer, number of lags for Newey-West correction
%
%   Output:
%       Q   - scalar, Newey-West standard error for x
%

    T = length(eps);           % Number of observations
    X = [ones(T,1), x];        % Design matrix with intercept

    % Diagonal matrix of squared residuals
    D = diag(eps.^2);
    % "A" matrix for NW formula (contemporaneous)
    A = (T-2) * (X' * D * X) / T;

    if L == 0
        % No autocorrelation adjustment
        Q = inv(X' * X) * A * inv(X' * X);
        Q = sqrt(Q(2,2));  % Return standard error for x (not intercept)
    else
        % Newey-West correction for autocorrelation up to lag L
        B = 0;
        for l = 1:L
            % Off-diagonal for lag l
            Dl = diag(eps((l+1):T) .* eps(1:(T-l)));
            wl = 1 - l/(L+1);   % Bartlett kernel weight
            % Add both "forward" and "backward" components
            B = B + wl * (X((l+1):T, :)' * Dl * X(1:(T-l), :) + ...
                         (X((l+1):T, :)' * Dl * X(1:(T-l), :))') / T;
        end
        B = (T-2) * B;

        Q = inv(X' * X) * (A + B) * inv(X' * X);
        Q = sqrt(Q(2,2)); % Return standard error for x
    end

end
