function [P0, P1, P1min] = NullAndAlt_var_bic(b, bias, se, s, alpha, bic)
%For lag length selection example Monte Carlo experiments, 
% this function returns:

%P0 - the null distribution of p-values under no p-hacking
%P1 - the distribution of p-values under thresholding p-hacking approach
%P1min - the distribution of p-values under minimum p-hacking approach

%The following inputs need to be provided:

%b - the vector of coefficients for corresponding regrssions
%bias - the vector of biases
%se - the vector of standard errors
%s = 1 or 2 for 1- or 2-sided tests respectively
%alpha - significance level used by researchers
%bic - the vector of BIC selected lags

M = length(bias(1, :));

P0 = zeros(M,1);
P1 = zeros(M,1);

t = (b + bias)./se;
if (s==1)
    p = 1 - normcdf(t);
end
if (s==2)
    p = 2*(1 - normcdf(abs(t)));
end



% if (up == 0)
% p = p(end:-1:1, :);
% end
for m = 1:M
P0(m,1) = p(bic(m), m)';
end

P1min = min(p)';

for m = 1:M
res = p(bic(m):end, m);

if (min(res)>alpha)
    P1(m) = min(res);
end
if (min(res)<=alpha)
    P1(m) = res(min(find(res<=alpha)));
end

end


end
