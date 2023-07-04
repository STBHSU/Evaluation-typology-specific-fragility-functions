function [im_out,success, theta_out, beta_out] = rtim_linear(k1, k0, im_haz, lambda_t, m_theta, c_theta, m_beta, c_beta, beta_haz, beta_cap, tol, max_iter, print)
%CALCULATE_IM_RISK : design IM to achieve design risk level
%   INPUTS:
%       k1, k0  : coefficients for linear hazard curve
%       im_haz  : the design im from hazard assessment
%       m_theta, b_theta    : coefficients for theta vs. design im
%       m_beta, c_beta    : coefficients for dispersion vs. design im
%       beta_haz: dispersion associated with hazard curve
%       beta_cap: dispersion of the median of the fragility curve
%       tol     : tolerance as fraction of target e.g. 0.001
%       max_iter: maximum iterations to perform e.g. 15
%       print   : == 1 and the output of each step is printed to screen
%       

% maximum acceptable error in the estimation of the risk-targeted im
max_error = lambda_t * tol;

% calculate the risk-targeted im
% preallocating for speed
ims = zeros([1,max_iter]);
theta = zeros([1,max_iter]);
beta = zeros([1,max_iter]);
lambda = zeros([1,max_iter]);
error = zeros([1,max_iter]);

% calculate the theta-based on y-intercept
theta_min = theta_linear(m_theta, c_theta, 0);
beta_min = sqrt((m_beta .* 0 + c_beta).^2 + beta_haz.^2 + beta_cap.^2);
lambda_min = linear_mafe(k0, k1, theta_min, beta_min);
error_min = lambda_t - lambda_min;

if error_min > 0
    % the assumed inherent lateral capacity dominates earthquake loading
    im_out = 0;
    success = 2;
    theta_out = theta_min;
    beta_out = beta_min;

else
    % earthquake loading assumed to dominate the assumed lateral capacity
    % create the first two points
    % alter first point until error is less than zero. This ensures convergence
    % second guess smaller to always head toward root
    guess = 0.5*im_haz;
    theta(1) = theta_linear(m_theta, c_theta, guess);
    beta(1) = sqrt((m_beta .* guess + c_beta).^2 + beta_haz.^2 + beta_cap.^2);
    lambda(1) = linear_mafe(k0, k1, theta(1), beta(1));
    error(1) = lambda_t - lambda(1);
    
    while error(1) > 0
        guess = 0.5 * guess;
        theta(1) = theta_linear(m_theta, c_theta, guess);
        beta(1) = sqrt((m_beta .* guess + c_beta).^2 + beta_haz.^2 + beta_cap.^2);
        lambda(1) = linear_mafe(k0, k1, theta(1), beta(1));
        error(1) = lambda_t - lambda(1);
    
    end
    
    ims(1:2) = [guess, 0.8*guess]; 
    
    theta(1:2) = theta_linear(m_theta, c_theta, ims(1:2));
    beta(1:2) = sqrt((m_beta .* ims(1:2) + c_beta).^2 + beta_haz.^2 + beta_cap.^2);
    lambda(1:2) = linear_mafe(k0, k1, theta(1:2), beta(1:2));
    error(1:2) = lambda_t - lambda(1:2);
    
    % printing status for debugging
    if print == true
        fprintf(['ii = %d\n' ...
                 'error = %f\n' ...
                 'pga: %f, theta: %f, beta: %f, lambda: %f\n\n'], ...
                 1, error(1), ims(1), theta(1), beta(1), lambda(1))
    end
    
    % solving for the risk-targeted PGA
    for ii = 2:max_iter
        % check if the ii lambda is suitable
        error(ii) = lambda_t - lambda(ii);
    
        % print output if required
        if print == true
                fprintf(['ii = %d\n' ...
                         'error = %f\n' ...
                         'pga: %f, theta: %f, beta: %f, lambda: %f\n\n'], ...
                         ii, error(ii), ims(ii), theta(ii), beta(ii), lambda(ii))
        end
        
        if abs(error(ii)) <= max_error
            break    
        else
            % use newtons method to get the new x value
            ims(ii+1) = newton(ims(ii), error(ii), ims(ii-1), error(ii-1));
            theta(ii+1) = theta_linear(m_theta, c_theta, ims(ii+1));
            beta(ii+1) = sqrt((m_beta .* ims(ii+1) + c_beta).^2 + beta_haz.^2 + beta_cap.^2);
            lambda(ii+1) = linear_mafe(k0, k1, theta(ii+1), beta(ii+1));
        end
    end
    
    if(ii < max_iter)
        im_out = ims(ii);
        success = 1;
        theta_out = theta(ii);
        beta_out = beta(ii);
    else
        im_out = ims(end);
        success = 0;
        theta_out = theta(ii);
        beta_out = beta(ii);
    end
end

end