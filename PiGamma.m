function [LogLikelihood] = PiGamma(X, a, b)

l = zeros(size(X));
h = exp(linspace(-15, 3, 1000));
Xs = linspace(-10,10, 1000)';
Mxs = normpdf(Xs-h);
Mhs = repmat(gampdf(h, a, b), length(Xs), 1);
Is = trapz(h, Mxs.*Mhs, 2);
Is = trapz(Xs,Is);

Mx = normpdf(X-h);
Mh = repmat(gampdf(h, a, b), length(X), 1);
I = trapz(h, Mx.*Mh, 2)/Is;
% parfor j = 1:length(X)
%     %I = integral(@(h) (normpdf(X(j)-h)).*gampdf(h, a, b), 0, Inf);
%     I = (normpdf(X(j)-h)).*gampdf(h, a, b);
%     I = trapz(h, I);
%     l(j) = log(I);
% end
LogLikelihood = sum(log(I));
end

