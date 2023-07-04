% Title: Comparing Prob. of Col. in DE using EC8 NA:2018
% Author: Nicholas Clemett
% Date: 10.10.22

% Description:
% This script calculates the values for the maps presented in Figure 14 of
% the paper
%   - the Pc is calculated for structures located in germany
%   - different typology-specific fragility models are used
%   - a piecewise linear relationship between median collapse pga and the design pga
%   - the dispersion of the hazard curve is NOT considered

clear
close all
clc

tic
fprintf("Running...\n")

%% Input Parameters

% file names for save and import files go here
site_data_file = "site_hazard_data_18.mat";
lit_file = "rts_specific_piecewise.csv";
save_name = "data_out\fig14_pc";

% list of indentifiers of fragility curves to use:
% to be able to plot the results this file needs to be run twice. Once with
%  "rc-mrf-m-pw" and once with "rc-mrf-m-rto".

ids = ["rc-mrf-m-pw"]; % --> the piecewise model with inherent capacity
% ids = ["rc-mrf-m-rto"]; % ---> regression through the origin

% hazard parameters
Occ = [0.10; 0.02]'; % Probability of exceedance
t = 50; % time span considered for occurence probabilities
incl_beta_haz = false; % include the dispersion due to hazard curve
pga_min = 0.01; % the min. 475yr PGA [ms-2] considered in calcs (>=0.01)


%% Precalculations
% check and create data_out folder
if not(isfolder(pwd + "\data_out"))
    mkdir(pwd + "\data_out")
end

% initialise hazard data
data = load(site_data_file);
fields = fieldnames(data);
sd = data.(fields{1});
sd = sd(find(sd.pga_475_median >= pga_min), :);
lat = sd.lat;
lon = sd.lon;

% redefine the maxcoord for the reduced list of coordinates
maxcoord = length(lat);

% get k0 and k1 parameters for the valid coordinates
lambda = -log(1 - Occ) / t; % mean annual frequency of occurence
[k0, k1] = linear_haz_params([sd.pga_475_median, sd.pga_2475_median], lambda);

%% Calculating the collapse risk using NA:2018 PGAs

% import the fragility curve data from literature file
lit_params = readtable(lit_file, "ReadRowNames", true);

for ii = 1:1:length(ids)
    fc = lit_params([ids(ii)],:);
    
    theta_min = fc.m_theta * fc.pga_plat;
    thetas = fc.m_theta .* sd.pga_475_median + fc.c_theta;
    thetas(thetas < theta_min) = theta_min;

    % additional dispersion due to variation of hazard curve
    if incl_beta_haz == true
        beta_hazs = abs(log(sd.pga_475_84) - log(sd.pga_475_median));
    else
        beta_hazs = zeros(length(lat),1);
    end

    beta_caps = fc.beta_cap .* ones(length(lat),1);
    betas = fc.m_beta .* ones(length(lat),1) .* sd.pga_475_mean + fc.c_beta;

    beta_tots = ((betas.^2) + (beta_caps.^2) + (beta_hazs.^2)).^(0.5);
    
    % collapse
    mafc = linear_mafe(k0, k1, thetas, beta_tots);
    
    % convert to probability of collapse in t years
    pc = (1 - exp(-mafc*t)) * 100; % in percent

    % calculate the moment estimators for the log normal distribution
    pc_median = exp(sum(log(pc)) / length(pc));
    pc_beta = sqrt((sum(log(pc / pc_median).^2)) / (length(pc) - 1));

    % saving variables
    save(save_name + "_" + fc.id + "_NA-2018" + ".mat", ...
    "lat", "lon", "k1", "k0", "mafc", "pc", "pc_median", "pc_beta")
end

fprintf("Done!\n")
toc

% plot_pc_pdfs








