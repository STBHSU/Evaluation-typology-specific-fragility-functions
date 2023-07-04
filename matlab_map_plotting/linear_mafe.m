function [lambda_ls] = linear_mafe(k0, k1, theta, beta)
%LINEAR_MAFE mean annual frequency of exceedance of a limit state
%   calculates the mean annual frequency of exceedence of a limit state
%   assuming that the fragility curve is a lognormal cdf with a median,
%   theta, and a dispersion, beta. The hazard curve is assumed to be linear
%   in log-log space with intercept k0 and slope k1

lambda_ls = k0 .* (theta.^(-k1)) .* exp(0.5 .* (k1.^2) .* (beta.^2));
end