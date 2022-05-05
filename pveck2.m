% function to examine all combinations of the z's as covariates and return
% the vector of pvalues from regressing y on x

function [bias, se] = pveck2(eps,x,z,k)

% note  k<dim(x)<15

[n,K]=size(z);

v=(1:1:K)';
s1=nchoosek(v,k);
%p=zeros(length(s1),1);
bias=zeros(length(s1),1);
se=zeros(length(s1),1);
e1 = zeros(1, k+2);
e1(1)=1;
for i=1:length(s1)
    
    x1=[x ones(n,1) z(:,s1(i,:))];
    %[b1,se]=ols1(y,x1,0,0);
    bias(i,1) = e1*inv(x1'*x1)*x1'*eps;
    eps_hat = x1*inv(x1'*x1)*x1'*eps - eps;
    se(i,1) = sqrt(e1*(inv(x1'*x1)*sum((eps_hat).^2)/(size(x1,1) - size(x1,2)))*e1');
    %t=b1(1,1)/sqrt(se(1,1));
    %p(i,1)=1-normcdf(t);    % upper tail rejections
end
