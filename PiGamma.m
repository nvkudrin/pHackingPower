function [LogLikelihood] = PiGamma(X, a, b)
%This function calculates the log-likelihood function for estimating
%parameters of Pi = Gamma(a,b) distribution
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
LogLikelihood = sum(log(I));
end

