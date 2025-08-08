% Program to generate Monte Carlo data for constructing p-hacked and non-p-hacked distributions.
% Paper: The Power of Tests for Detecting p-Hacking
% Authors: G. Elliott, N. Kudrin, K. Wuthrich

clear all
addpath('Functions')
rng(123)

nobs = 200;      % Sample size used by researchers
nmc  = 1e6;      % Number of draws for constructing distributions

lev = [1, 2, 4, 5, 10];   % Cluster levels for Example 3.2

% Preallocate storage for results
mcout_bias_cov     = []; % Bias terms covariate selection
mcout_se_cov       = []; % SE terms covariate selection
mcout_bias_iv5     = []; % Bias terms IV selection (5 IVs)
mcout_se_iv5       = []; % SE terms IV selection (5 IVs)
mcout_F_iv5        = []; % First stage F for IV selection (5 IVs)
mcout_bias_iv3     = []; % Bias terms IV selection (3 IVs)
mcout_se_iv3       = []; % SE terms IV selection (3 IVs)
mcout_F_iv3        = []; % First stage F for IV selection (3 IVs)
mcout_bias_laglength     = []; % Bias terms for lag length selection
mcout_se_laglength       = []; % SE terms for lag length selection
mcout_bic_lag      = []; % BIC-selected lag length
mcout_bias_cluster = []; % Bias terms for cluster level selection
mcout_se_cluster   = []; % SE terms for cluster level selection

sc = parallel.pool.Constant(RandStream('mlfg6331_64'));
tic

parfor imc = 1:nmc

    stream = sc.Value;           % Extract parallel stream
    stream.Substream = imc;

    %%%%%%% Covariate selection
    ks = 7;  % Maximum number of controls to select between

    x_cov   = randn(stream, nobs, 1);
    gamm = rand(stream, ks, 1) * 1.6 - 0.8;         % gamma ~ U[-0.8, 0.8]

    % Generate controls correlated with x1 and independent of each other
    zz = zeros(nobs, ks);        
    for iz = 1:ks
        u = randn(stream, nobs, 1) * sqrt(1 - gamm(iz)^2);  
        zz(:, iz) = gamm(iz) * x_cov + u;
    end

    eps_cov = randn(stream, nobs, 1);   % Regression error
    Bias_cov = [];
    StdErr_cov = [];
    for j = 0:ks 
        [B_cov, SE_cov] = RegSpecifications(eps_cov, x_cov, zz, ks - j);
        Bias_cov   = [Bias_cov; B_cov];
        StdErr_cov = [StdErr_cov; SE_cov];
    end
    mcout_bias_cov = [mcout_bias_cov, Bias_cov];      % Record bias terms
    mcout_se_cov   = [mcout_se_cov, StdErr_cov];      % Record standard errors

    %%%%%%% IV selection (K=5)
    IV5    = zz(:, 1:5);
    eps_iv = randn(stream, nobs, 1);
    x_iv5  = IV5 * (1 + 2 * rand(stream, 5, 1)) + 0.5 * eps_iv + sqrt(1 - 0.5^2) * randn(stream, nobs, 1);

    Bias_iv5   = [];
    StdErr_iv5 = [];
    Fstat5     = [];
    for j = 0:4
        [Biv, SEiv, Fstat_5] = IVSpecifications(eps_iv, x_iv5, IV5, 5 - j);
        Bias_iv5   = [Bias_iv5; Biv];
        StdErr_iv5 = [StdErr_iv5; SEiv];
        Fstat5     = [Fstat5; Fstat_5];
    end
    mcout_bias_iv5 = [mcout_bias_iv5, Bias_iv5];
    mcout_se_iv5   = [mcout_se_iv5, StdErr_iv5];
    mcout_F_iv5    = [mcout_F_iv5, Fstat5];

    %%%%%%% IV selection (K=3)
    IV3    = zz(:, 1:3);
    x_iv3  = IV3 * (1 + 2 * rand(stream, 3, 1)) + 0.5 * eps_iv + sqrt(1 - 0.5^2) * randn(stream, nobs, 1);

    Bias_iv3   = [];
    StdErr_iv3 = [];
    Fstat3     = [];
    for j = 0:2
        [Biv3, SEiv3, Fstat_3] = IVSpecifications(eps_iv, x_iv3, IV3, 3 - j);
        Bias_iv3   = [Bias_iv3; Biv3];
        StdErr_iv3 = [StdErr_iv3; SEiv3];
        Fstat3     = [Fstat3; Fstat_3];
    end
    mcout_bias_iv3 = [mcout_bias_iv3, Bias_iv3];
    mcout_se_iv3   = [mcout_se_iv3, StdErr_iv3];
    mcout_F_iv3    = [mcout_F_iv3, Fstat3];

    %%%%%%% Variance manipulation: lag length selection
    x_var   = randn(stream, nobs, 1);
    eps_var = randn(stream, nobs, 1);
    X_var   = [x_var ones(nobs, 1)];
    e1 = zeros(1, 2); e1(1) = 1;

    bias_var_full = (X_var' * X_var) \ (X_var' * eps_var); % Full bias vector
    bias_var = e1 * bias_var_full;                         % Bias term for x_var

    e_var  = eps_var - X_var * bias_var_full;
    eX_var = e_var .* x_var;                               % Input for BIC

    mcout_bic_lag   = [mcout_bic_lag, LagsBIC(eX_var)];                 % Record BICs
    mcout_bias_laglength  = [mcout_bias_laglength, bias_var * ones(5, 1)];          % Record bias terms

    % Calculate Newey-West SE for different lags
    se_var = NeweyWest(e_var, x_var, 0);
    for L = 1:4
        se_var = [se_var; NeweyWest(e_var, x_var, L)];
    end
    mcout_se_laglength = [mcout_se_laglength, se_var];

    %%%%%%% Variance manipulation: clustering level
    x_clust   = randn(stream, nobs, 1);
    eps_clust = randn(stream, nobs, 1);
    X_clust   = [x_clust ones(nobs, 1)];
    e1 = zeros(1, 2); e1(1) = 1;
    bias_clust = e1 * ((X_clust' * X_clust) \ (X_clust' * eps_clust));
    mcout_bias_cluster = [mcout_bias_cluster, bias_clust * ones(length(lev), 1)];

    eps_hat_clust = X_clust * ((X_clust' * X_clust) \ (X_clust' * eps_clust)) - eps_clust;
    se_clust = [];
    for g = 1:length(lev)
        G   = nobs / lev(g); % number of clusters
        u_g = reshape(eps_hat_clust, lev(g), G);
        X_g = reshape(X_clust', 2, lev(g), G);
        Ome = 0;
        for j = 1:G
            Ome = Ome + X_g(:, :, j) * u_g(:, j) * u_g(:, j)' * X_g(:, :, j)';
        end
        Ome = G * (nobs - 1) * inv(x_clust' * x_clust) * Ome * inv(x_clust' * x_clust) / (G - 1) / (nobs - 2);
        se_clust = [se_clust; sqrt(Ome(1, 1))];
    end
    mcout_se_cluster = [mcout_se_cluster, se_clust];
end

toc
s = rng;
mkdir('Raw_Pseudo_Data')
save('Raw_Pseudo_Data/MC_raw_pseudo_data.mat')
