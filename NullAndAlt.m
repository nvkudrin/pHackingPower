function [P0, P1,P1min,  Bias0, Bias1, Bias1min] = NullAndAlt(b, bias, se, s, alpha, Ktotal, K)
%For covariate selection and IV selection Monte Carlo experiments, 
% this function returns:

%P0 - the null distribution of p-values under no p-hacking
%P1 - the distribution of p-values under thresholding p-hacking approach
%P1min - the distribution of p-values under minimum p-hacking approach
%Bias0 - the average bias under no p-hacking
%Bias1 - the average bias under thresholding p-hacking approach
%Bias1min - the average bias under minimum p-hacking approach

%The following inputs need to be provided:

%b - the vector of coefficients for corresponding regrssions
%bias - the vector of biases
%se - the vector of standard errors
%s = 1 or 2 for 1- or 2-sided tests respectively
%alpha - significance level used by researchers
%Ktotal = total number of controls/instruments used in simulation
%K - number of controls/instruments researchers have access to

M = length(bias(1, :));
P0 = zeros(M,1);
P1 = zeros(M,1);
Bias0 = bias(1, :)';
Bias1 = zeros(M,1);
Bias1min = zeros(M,1);

v = (1:1:Ktotal);
V = [];
for j = 0:(Ktotal-1)
    a = nchoosek(v,Ktotal-j);
    V = [V;(max(a, [], 2)<=K)];
end
V = logical(V);
b = b(V, :);
bias = bias(V, :);
se = se(V, :);
t = (b + bias)./se;
if (s==1)
    p = 1 - normcdf(t);
end
if (s==2)
    p = 2*(1 - normcdf(abs(t)));
end

P0 = p(1, :)';
%P1min = min(p)';

[P1min, I] = min(p, [], 1);
P1min = P1min';



Ind = 2.^[0:K];
for m = 1:M
Bias1min(m, 1) = bias(I(m), m);    
X = p(:, m);
res = zeros(1,K);
for j = 1:(K)
    res(j) = min(X(Ind(j): (Ind(j+1)-1)));
end
if (min(res)>alpha)
    P1(m) = min(res);
    Bias1(m) = bias(min(find(X == P1(m))), m);
end
if (min(res)<=alpha)
    P1(m) = res(min(find(res<=alpha)));
    Bias1(m) = bias(min(find(X == P1(m))), m);
end

end


end

