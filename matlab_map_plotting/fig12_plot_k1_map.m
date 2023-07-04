clear
close all
clc

%% INPUTS

% file containing variables for mapping
save_file = "fig12_k1-parameter-map";

figure_folder = pwd;

% figure size, aspect ratio and offsets
fig_h = 8;
img_ar = 0.83; 
lft = 1.5; 
gap = 2.0;
rgt = 2.0;  
bot = 1.0;  
top = 0.2; 

% Font Settings
font = "Times";
label_fs = 12;
tick_fs = 9;

% Import hazard data for plotting
pga_data = readmatrix("pga.csv");
lat = pga_data(:,1);
lon = pga_data(:,2);
% pga_median = pga_data(:,[3,6,9]);
pga_median = pga_data(:,[3,9]);

bad_coords = find(pga_median(:,1) == 0.0);
f_lat = lat;
f_lat(bad_coords) = [];
f_lon = lon;
f_lon(bad_coords) = [];
f_pga_median = pga_median;
f_pga_median(bad_coords, :) = [];

% Auftretenswahrscheinlichkeit - This is converted to MAFE for calculations
% Occ = [0.10; 0.05; 0.02]';
Occ = [0.10; 0.02]';
t = 50; % time span considered for occurence probabilities
lambda = -log(1 - Occ) / t; % Mean annual frequency of occurence

% Calculate the k0 and k1 parameters
[k0, k1] = linear_haz_params(f_pga_median, lambda);


%% PLOTTING

% create 2D matrix of points for contour plot
maxcoord = length(f_lat);
delt_lat = 0.1;
delt_lon = 0.1;
y = min(f_lat):delt_lat:max(f_lat);
x = min(f_lon):delt_lon:max(f_lon);
map_k0 = zeros(max(size(y)),max(size(x)));
map_k1 = zeros(max(size(y)),max(size(x)));

k0_max = max(k0);
k0_min = min(k0);
k1_max = max(k1);
k1_min = min(k1);

for coord=1:maxcoord
        p = 1+ round((f_lon(coord)-min(f_lon))/delt_lon);
        q = 1+ round((f_lat(coord)-min(f_lat))/delt_lat);
        
        map_k0(q,p) = k0(coord);
        map_k1 (q,p) = k1(coord);
end

map_k0(map_k0 == 0) = NaN;
map_k1(map_k1==0)=NaN;

% Create a new figure for the maps
% image and figure sizes
img_h = fig_h - (top + bot);
img_w = img_h * img_ar;
fig_w = img_w + lft + rgt;

% position vectors for two images side by side
% pos1 = [lft, bot, img_w, img_h];
% pos2 = [lft + img_w + gap, bot, img_w, img_h];
pos2 = [lft, bot, img_w, img_h];

% import border points
border = readmatrix("new_border.csv");

% Create Figure
f=figure("Units","centimeters", "Position",[10,10,fig_w,fig_h]);

% p1 = subplot(1,2,1); % first subplot for risk-based pga
% hold on
% contourf(x, y, map_k0)%,[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])
% plot(border(:,1), border(:,2), "LineWidth",2, "Color","black")
% hold off
% cb = colorbar("FontName", font, "Fontsize", tick_fs);
% cb.Label.String = "k0-Values"; 
% colormap(p1, "jet")
% %caxis([0, 2.0])
% ax = gca;
% set(ax,'Units','centimeters', "Position", pos1)
% ax.XAxis.FontSize = tick_fs;
% ax.XAxis.FontName = font;
% ax.YAxis.FontSize = tick_fs;
% ax.YAxis.FontName = font;
% set(get(p1, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
% set(get(p1, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

p2 = subplot(1,2,2); % second subplot is the risk-coefficient
hold on
contourf(x, y, map_k1, [0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2])
plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
hold off
caxis([0.8, 2.2])
cb = colorbar("FontName", font, "Fontsize", tick_fs);
cb.Label.String = "k1-Values";
colormap(p2, "jet")
ax = gca;
ax.XAxis.FontSize = tick_fs;
ax.XAxis.FontName = font;
ax.YAxis.FontSize = tick_fs;
ax.YAxis.FontName = font;
set(ax,'Units','centimeters','Position', pos2)
set(get(p2, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
set(get(p2, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

annotation(textbox=[0.36, 0.56, 0.3, 0.3], ...
           String=("Max. k_{1} = " + string(round(k1_max,3)) + newline + ...
                   "Min. k_{1} = " + string(round(k1_min,3))), ...
           FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);

saveas(f, figure_folder + "\" + save_file + ".png")
saveas(f, figure_folder + "\" + save_file + ".pdf")
