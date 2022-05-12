% This code implements Newey-West standard errors used in the lag length
% selection example
% Authors: G. Elliott, N. Kudrin, K. Wuthrich
function [Q] = NeweyWest(eps,x,L)
% Inputs: eps - model error term, x - regressor, L - the number of lags
% Output: Q - NW standard error
T = length(eps);
X = [ones(T,1), x];
D = diag(eps.^2);
A = (T-2)*(X'*D*X)/T;
if (L==0)
    Q = inv(X'*X)*A*inv(X'*X);
    Q = sqrt(Q(2,2));
end
if (L>0)
B = 0;
for l = 1:L
   Dl = diag(eps((l+1):T).*eps(1:(T-l)));
   wl = 1 - l/(L+1);
   B = B + wl*((X((l+1):T, :)'*Dl*X(1:(T-l), :)) + (X((l+1):T, :)'*Dl*X(1:(T-l), :))')/T;
end
B = (T-2)*B;

Q = inv(X'*X)*(A+B)*inv(X'*X);
Q = sqrt(Q(2,2));
end