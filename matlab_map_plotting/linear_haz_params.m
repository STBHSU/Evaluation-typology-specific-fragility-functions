function [k0,k1] = linear_haz_params(ims,lambdas, varargin)
%LINEAR_HAZ_PARAMS: calculate the k0 and k1 parameters for hazard curves
%   This function calculates the k0 and k1 parameters for the hazard curve
%   based on a linear approximation in log-log space.
%   The linear fit is determine through least sqaures regression of the
%   given points.
%   ASSUMES: that all of the ims are non-zero
%
%INPUT ARGUMENTS:
%   ims     :   an mxn matrix of the pga values. m different sites, n 
%               points from the hazard curve. It is assumed that the
%               columns are order left to right by reducing MAFE
%   lambdas  :  a 1xn array of the MAFE for each column in ims
%   lin_type :  optional. String. The type of linear fit. 
%               - "regression"
%               - "two-point"
%               - "jalayer03"
%
%OUTPUTS:
%   k0      :   the y intercept of the linear fit
%   k1      :   the slope of the linear fit


maxcoord = length(ims(:,1));
k0 = zeros([maxcoord, 1]);
k1 = zeros([maxcoord, 1]);

for coord = 1:maxcoord
    a = log(lambdas);
    b = log(ims(coord,:));
    f = polyfit(b,a,1);         % f(1)=Steigung,  f(2)=Achsenabschnitt
    k0(coord) = exp(f(2));
    k1(coord)  = f(1) * -1;
end

end
    