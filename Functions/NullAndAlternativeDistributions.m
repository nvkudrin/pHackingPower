function [P0, P1, P1min, Bias0, Bias1, Bias1min] = NullAndAlternativeDistributions(b, bias, se, s, alpha, Ktotal, K, GeneralToSpecific, example, Fstat)
%NULLANDALTERNATIVEDISTRIBUTIONS
%   For covariate selection and IV selection Monte Carlo experiments, returns:
%     P0      - null distribution of p-values (no p-hacking)
%     P1      - p-value distribution under thresholding p-hacking
%     P1min   - p-value distribution under minimum p-hacking
%     Bias0   - bias under no p-hacking
%     Bias1   - bias under thresholding p-hacking
%     Bias1min- bias under minimum p-hacking
%
%   Inputs:
%     b                - [matrix] regression coefficients
%     bias             - [matrix] bias values
%     se               - [matrix] standard errors
%     s                - 1 or 2 (one- or two-sided test)
%     alpha            - significance level
%     Ktotal           - total number of controls/instruments in simulation
%     K                - number of controls/instruments researchers use
%     GeneralToSpecific- 1 for general-to-specific, 0 otherwise
%     example          - 'CovSel' for Covariate Selection, 'IV' for IV
%     Fstat            - [] for no first-stage screening, or [matrix] first-stage F stats for IV
%
%   Outputs:
%     P0, P1, P1min    - [M x 1] p-value distributions for M MC draws
%     Bias0, Bias1, Bias1min - [M x 1] bias values for each scenario

    M = size(bias, 2);           % Number of Monte Carlo draws

    P1     = zeros(M, 1);        % p-hacked distribution (thresholding)
    Bias0  = bias(1, :)';        % Bias under no p-hacking
    Bias1  = zeros(M, 1);        % Bias under thresholding p-hacking
    Bias1min = zeros(M, 1);      % Bias under minimum p-hacking

    % Identify all relevant specifications for the given K
    v = 1:Ktotal;
    V = [];
    for j = 0:(Ktotal-1)
        a = nchoosek(v, Ktotal-j);
        V = [V; (max(a, [], 2) <= K)];
    end

    if example == "CovSel"
        V = logical([V; 1]);
    elseif example == "IV"
        V = logical(V);
    end

    % Keep only relevant specs
    b    = b(V, :);
    bias = bias(V, :);
    se   = se(V, :);

    % Compute t-stats and p-values
    t = (b + bias) ./ se;
    if s == 1
        p = 1 - normcdf(t);                  % One-sided
    elseif s == 2
        p = 2 * (1 - normcdf(abs(t)));       % Two-sided
    else
        error('Input s must be 1 or 2');
    end

    % No p-hacking: first/last spec depending on direction
    if GeneralToSpecific == 1
        P0 = p(1, :)';
    else
        P0 = p(end, :)';
    end

    % If Fstat provided, exclude p-values from weak first stages (IV only)
    if ~isempty(Fstat)
        for m = 1:M
            if max(Fstat(:, m)) > 10
                p(Fstat(:, m) < 10, m) = 10; % Set to 10 (never selected)
            end
        end
    end

    % Minimum p-hacking
    [P1min, I] = min(p, [], 1);
    P1min = P1min';
    
    % Record bias under minimum p-hacking
    for m = 1:M
        Bias1min(m, 1) = bias(I(m), m);
    end

    % Identify p-hacking step sizes for thresholding
    Ind = [];
    for j = 0:(K-1)
        Ind = [Ind, nchoosek(K, K-j)];
    end
    if example == "CovSel"
        Ind = [Ind, 1];
    end
    Ind = cumsum(Ind);
    Ind_size = length(Ind);

    % Thresholding p-hacking (stepwise min-by-step)
    for m = 1:M
        X = p(:, m);
        % Direction of search
        if GeneralToSpecific == 0
            X = flip(X);
        end
        res = zeros(1, Ind_size);
        res(1) = X(1);
        for j = 2:Ind_size
            res(j) = min(X(Ind(j-1)+1 : Ind(j)));
        end

        % Select p-value and bias as per thresholding strategy
        if min(res) > alpha
            P1(m) = min(res);
            if GeneralToSpecific == 1
                Bias1(m) = bias(min(find(X == P1(m))), m);
            else
                Bias1(m) = bias(min(find(flip(X) == P1(m))), m);
            end
        else
            idx = min(find(res <= alpha));
            P1(m) = res(idx);
            if GeneralToSpecific == 1
                Bias1(m) = bias(min(find(X == P1(m))), m);
            else
                Bias1(m) = bias(min(find(flip(X) == P1(m))), m);
            end
        end
    end

end
