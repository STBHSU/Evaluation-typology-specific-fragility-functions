% Title: Risk-targeted PGA using generic fragility curves
% Author: Nicholas Clemett
% Date: 07.10.22

% Description:
%   - The risk-targeted PGA is calculated using the linear hazard model
%   - different fragility models are used from the literature
%   - a linear relationship between median collapse pga and the design pga
%   - the dispersion of the hazard curve is not considered

clear
close all
clc

%% Input Parameters
tic
fprintf("Running...\n")

% filenames
site_data_file = "site_hazard_data_18.mat";
lit_file = "rts_specific.csv";
save_name = pwd;
figure_name = "rt_pga_map_median";
fig_folder = pwd;

% list of indentifiers of fragility curves to use
ids = ["rc-mrf-m-rto"];
% ids = ["rc-mrf-m"];
% ids = ["rc-mrf-r" "s-mrf" "rc-wds"];

% Auftretenswahrscheinlichkeit - This is converted to MAFE for calculations
Occ = [0.10; 0.02]';
t = 50; % time span considered for occurence probabilities

% parameters for determining risk-targeted PGA
lambda_t = 0.00005;  % target collapse probability: 0.01 in 50 years
tol = 0.005;   % convergence tolerance
max_iter = 50; % seems reasonable
print = false;  % print output of interations

% parameters for determining the appropriate fragility curves
% c_theta = 0;  % intercept for theta equation
% m_beta = 0;  % slope for the dispersion equation, =0 for constant
% beta_cap = 0.0; % additional dispersion due to variation of median value
incl_beta_haz = false; % considered increased dispersion from hazard curve

%other
pga_min = 0.01; % the min. 475yr PGA [ms-2] used in calcs (>=0.01)

% import the site hazard data and filter based on minimum PGA
data = load(site_data_file);
fields = fieldnames(data);
sd = data.(fields{1});
sd = sd(find(sd.pga_475_median >= pga_min), :);

% import the fragility curve data from literature file
lit_params = readtable(lit_file, "ReadRowNames", true);

%% Precalculations
maxcoord = length(sd.lat);

% Mean annual frequency of occurence
lambda = -log(1 - Occ) / t;

% preallocations
pga_risk = zeros(maxcoord,1);
beta_tot = zeros(maxcoord,1);
converged = zeros(maxcoord,1);
unconverged = zeros(maxcoord,1);
too_strong = zeros(maxcoord,1);

% get k0 and k1 parameters
[sd.k0, sd.k1] = linear_haz_params([sd.pga_475_median, sd.pga_2475_median], lambda);

%% Calculating the risk-targeted PGA

% looping for each fragility curve selected
for ii = 1:1:length(ids)
    fc = lit_params([ids(ii)],:);
    
%     m_theta = exp(-norminv(fc.pcx) * fc.beta);
%     c_beta = fc.beta;
    
    converged_count = 0;
    unconverged_count = 0;
    too_strong_count = 0;
    
    fprintf("Calculating Map Values...")
    
    for coord = 1:maxcoord

        s = sd(coord, :);
%         fprintf("lat: %f, lon: %f\n", s.lat, s.lon)
        
        % additional dispersion due to variation of hazard curve
        if incl_beta_haz == true
            beta_haz = abs(log(s.pga_475_84) - log(s.pga_475_median));
        else
            beta_haz = 0;
        end
    
        % find the risk-targeted design pga for one site
        [pga_risk(coord), scs, ~, beta_tot(coord)] = rtim_linear(s.k1, s.k0, s.pga_475_median, ...
                                                  lambda_t, fc.m_theta, fc.c_theta, fc.m_beta, fc.c_beta, ...
                                                  beta_haz, fc.beta_cap, tol, max_iter, print); 
    
        if scs == 1
            %fprintf('coordinate %d converged successfully\n', coord);
            converged_count = converged_count + 1;
            converged(converged_count) = coord;
        elseif scs == 0
            %fprintf('coordinate %d DID NOT converge\n', coord);
            unconverged_count = unconverged_count + 1;
            unconverged(unconverged_count) = coord;
        elseif scs == 2
            too_strong_count = too_strong_count + 1;
            too_strong(too_strong_count) = coord;
        end
    
    end
    fprintf("Done!\n")

    % trimming zeros
    converged = converged(converged ~= 0);
    unconverged = unconverged(unconverged ~= 0);
    too_strong = too_strong(too_strong ~= 0);
    
    if unconverged_count > 0
        fprintf("WARNING: %d coordinates did not converge! - Check outputs\n", ...
                unconverged_count)
    else
        fprintf("All coordinates converged successfully\n")
    end
    
    % Risk Koeffizient Cr
    Cr = pga_risk ./ sd.pga_475_median;

    % saving variables
    k1 = sd.k1; k0 = sd.k0; lat = sd.lat; lon = sd.lon;
    file_name = save_name + "_" + fc.id + ".mat";
    save(file_name, "lat", "lon", "Cr", "pga_risk", "k1", "k0")

    plot_rt_pga_maps(figure_name + "_" + fc.id, fig_folder, lat, lon, Cr, pga_risk, 0.1, 0.1, ...
                     fc.id + ": Risk-targeted PGA and Risk Coefficient")

end








