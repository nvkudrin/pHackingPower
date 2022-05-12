% This code simulates data for Monte Carlo analysis from "(When) Can we detect p-hacking"
% Authors: G. Elliott, N. Kudrin, K. Wuthrich
%%
clear all
rng(123)
nobs=200; %sample sizes
nmc = 1000000; %the number of MC replications

% Variables to store the results
mcout_bias = [];
mcout_se = [];
mcout_bias_iv5 = [];
mcout_se_iv5 = [];
mcout_F_iv5 = [];
mcout_bias_iv3 = [];
mcout_se_iv3 = [];
mcout_F_iv3 = [];
mcout_bias_var = [];
mcout_se_var = [];
mcout_bic_var = [];

sc = parallel.pool.Constant(RandStream('mlfg6331_64'));
tic;
jj=0;
parfor imc=1:nmc

    stream = sc.Value;        % Extract the stream from the Constant
    stream.Substream = imc;
    % EXAMPLE 1: Covariate selection

% x1 is variable of interest
% z's are additional controls, correlated with x1 and indep of each other

ks=7;  % maximum number of regressors to select between

x1 = randn(stream, nobs,1);
gamm=rand(stream, ks,1)*1.6-0.8;         % generate gamm uniform on [-0.8,0.8]


zz=zeros(nobs,ks);        % preallocate space for z
for iz=1:ks
    u=randn(stream, nobs,1)*sqrt(1-gamm(iz,1)*gamm(iz,1));  
    zz(:,iz)=gamm(iz,1)*x1+u;
end

eps = randn(stream, nobs,1); % error term


    X1=[x1 ones(nobs,1) zz];

 
    e1 = zeros(1, 7+2);
    e1(1)=1;
    bias = e1*inv(X1'*X1)*X1'*eps;
    eps_hat = X1*inv(X1'*X1)*X1'*eps - eps;
    s7 = sqrt(e1*(inv(X1'*X1)*sum((eps_hat).^2)/(nobs-9))*e1');

    Bias = bias;
    StdErr = s7;
    for j = 1:(ks-1) 
        [B, SE] = pveck2(eps, x1, zz, ks-j);
    Bias = [Bias;B];
    StdErr = [StdErr; SE];
    end
    mcout_bias = [mcout_bias,Bias];
    mcout_se = [mcout_se,StdErr];




% EXAMPLE 2: IV selection (K=5)
IV = zz(:,1:5);

eps_iv = randn(stream, nobs, 1);
x_iv = IV*(1+2*rand(stream, 5, 1)) + 0.5*eps_iv+sqrt(1-0.5^2)*randn(stream, nobs, 1);
x5 = [ones(nobs,1) x_iv];
z_iv = [ones(nobs, 1) IV];
Pz = z_iv*inv(z_iv'*z_iv)*z_iv';

pihat5 = inv(z_iv'*z_iv)*z_iv'*x_iv;
V5 = (x_iv'*(eye(nobs) - Pz)*x_iv/(nobs - 6))*inv(z_iv'*z_iv);
Fstat5 = pihat5(2:end)'*inv(V5(2:end, 2:end))*pihat5(2:end)/5;
mcout_F_iv5 = [mcout_F_iv5, Fstat5];

    e1 = zeros(1, 2);
    e1(2)=1;
    bias_iv = e1*inv(x5'*Pz*x5)*(x5'*Pz*eps_iv);
    eps_hat_iv = x5*inv(x5'*Pz*x5)*x5'*Pz*eps_iv - eps_iv;
    s5 = sqrt(e1*(inv(x5'*Pz*x5)*sum((eps_hat_iv).^2)/(nobs-2))*e1'); %divide by n?

    Bias_iv = bias_iv;
    StdErr_iv = s5;
    for j = 1:(5-1) 
        [Biv, SEiv] = pveckIV(eps_iv, x_iv, IV, 5-j);
    Bias_iv = [Bias_iv;Biv];
    StdErr_iv = [StdErr_iv; SEiv];
    end
    mcout_bias_iv5 = [mcout_bias_iv5,Bias_iv];
    mcout_se_iv5 = [mcout_se_iv5,StdErr_iv];
    
    % EXAMPLE 2: IV selection (K=5)
IV3 = zz(:,1:3);

eps_iv3 = randn(stream, nobs, 1);
x_iv3 = IV3*(1+2*rand(stream, 3, 1)) + 0.5*eps_iv3+sqrt(1-0.5^2)*randn(stream, nobs, 1);
x3 = [ones(nobs,1) x_iv3];
z_iv3 = [ones(nobs, 1) IV3];
Pz3 = z_iv3*inv(z_iv3'*z_iv3)*z_iv3';

pihat3 = inv(z_iv3'*z_iv3)*z_iv3'*x_iv3;
V3 = (x_iv3'*(eye(nobs) - Pz3)*x_iv3/(nobs - 4))*inv(z_iv3'*z_iv3);
Fstat3 = pihat3(2:end)'*inv(V3(2:end, 2:end))*pihat3(2:end)/3;
mcout_F_iv3 = [mcout_F_iv3, Fstat3];

    e13 = zeros(1, 2);
    e13(2)=1;
    bias_iv3 = e13*inv(x3'*Pz3*x3)*(x3'*Pz3*eps_iv3);
    eps_hat_iv3 = x3*inv(x3'*Pz3*x3)*x3'*Pz3*eps_iv3 - eps_iv3;
    s3 = sqrt(e13*(inv(x3'*Pz3*x3)*sum((eps_hat_iv3).^2)/(nobs-2))*e13'); %divide by n?

    Bias_iv3 = bias_iv3;
    StdErr_iv3 = s3;
    for j = 1:(3-1) 
        [Biv3, SEiv3] = pveckIV(eps_iv3, x_iv3, IV3, 3-j);
    Bias_iv3 = [Bias_iv3;Biv3];
    StdErr_iv3 = [StdErr_iv3; SEiv3];
    end
    mcout_bias_iv3 = [mcout_bias_iv3,Bias_iv3];
    mcout_se_iv3 = [mcout_se_iv3,StdErr_iv3];
    
    
    
    % EXAMPLE 3: Variance manipulation
    
    x_var = randn(stream, nobs, 1);
    eps_var = randn(stream, nobs, 1);
    x_var2=[x_var ones(nobs,1)];
    e1 = zeros(1, 2);
    e1(1)=1;
    bias_var = e1*inv(x_var2'*x_var2)*x_var2'*eps_var;
    
    eX_var = (x_var*bias_var+eps_var).*x_var;
    mcout_bic_var = [mcout_bic_var, BIC(eX_var)];

    mcout_bias_var = [mcout_bias_var, bias_var*ones(5, 1)];
    
    se_var = NeweyWest(eps_var, x_var, 0);
    for L = 1:4
    se_var = [se_var; NeweyWest(eps_var, x_var, L)];
    end
    mcout_se_var = [mcout_se_var, se_var];
    
end
toc;
 s = rng;
 save('DGPs/MC_distributions_Apr6.mat')
%end
