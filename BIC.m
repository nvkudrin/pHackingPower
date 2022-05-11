% This code selects the initial number of lags for the lag length selection example
% using BIC criterion
% Authors: G. Elliott, N. Kudrin, K. Wuthrich
%%
function [k] = BIC(y)
K = [1,2,3,4,5];
T = length(y);
bic = zeros(1,5);
for j = 1:5
    Y = [];
    for i = 1:j
  Y = [Y, y((j+1-i):(T-i+1))];  
    end
X = [ones(length(Y(:,1)), 1), Y(:, 2:end)];
M = eye(length(Y(:,1))) - X*inv(X'*X)*X';
e = M*Y(:,1); 
t = length(e);
bic(j) = log(sum(e.^2)/t)+(j+1)*log(T)/T;%+t + t*log(2*pi);
end
k = find(bic == min(bic));
end

