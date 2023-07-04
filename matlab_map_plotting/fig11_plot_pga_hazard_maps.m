clear
close all
clc

%% INPUTS

% file containing variables for mapping
save_file = "Fig11_pga-hazard-map";

figure_folder = pwd;

% figure size, aspect ratio and offsets
fig_w = 16;
img_ar = 0.83; 
lft = 1.5; 
gap = 1.0;
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
pga_475_median = pga_data(:,3);
pga_2475_median = pga_data(:,9);

%% PLOTTING

% create 2D matrix of points for contour plot
maxcoord = length(lat);
delt_lat =0.1;
delt_lon =0.1;
y = min(lat):delt_lat:max(lat);
x = min(lon):delt_lon:max(lon);
map_475 = zeros(max(size(y)),max(size(x)));
map_2475 = zeros(max(size(y)),max(size(x)));

PGA_475_median_max = max(pga_475_median);
PGA_475_median_min = min(pga_475_median);
PGA_2475_median_max = max(pga_2475_median);
PGA_2475_median_min = min(pga_2475_median);

for coord=1:maxcoord
        p = 1+ round((lon(coord)-min(lon))/delt_lon);
        q = 1+ round((lat(coord)-min(lat))/delt_lat);
        
        map_475(q,p) = pga_475_median(coord);
        map_2475 (q,p) = pga_2475_median(coord);
end

map_475(map_475 == 0) = NaN;
map_2475(map_2475==0)=NaN;

% Create a new figure for the maps
% image and figure sizes
img_w = (fig_w - (lft + gap + rgt)) / 2;
img_h = img_w / img_ar;
fig_h = img_h + bot + top;

% position vectors for two images side by side
pos1 = [lft, bot, img_w, img_h];
pos2 = [lft + img_w + gap, bot, img_w, img_h];

% import border points
border = readmatrix("new_border.csv");

% Create Figure
f=figure("Units","centimeters", "Position",[10,10,fig_w,fig_h]);

p1 = subplot(1,2,1); % first subplot 475 RTP PGA
hold on
contourf(x, y, map_475,[0:0.2:3.2])
plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
hold off
% colorbar("FontName", font, "Fontsize", tick_fs) 
colormap(p1, "jet")
caxis([0, 3.0])
ax = gca;
set(ax,'Units','centimeters', "Position", pos1)
ax.XAxis.FontSize = tick_fs;
ax.XAxis.FontName = font;
ax.YAxis.FontSize = tick_fs;
ax.YAxis.FontName = font;
set(get(p1, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
set(get(p1, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

annotation(textbox=[0.18, 0.58, 0.3, 0.3], ...
           String=("Max. a_{g,475} = " + string(round(PGA_475_median_max,3)) + " m·s^{-2}" + newline + ...
                   "Min. a_{g,475} = " + string(round(PGA_475_median_min,3)) + " m·s^{-2}"), ...
           FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);

annotation(textbox=[0.09, 0.88, 0.1, 0.1], String="(a)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

p2 = subplot(1,2,2); % second subplot 2475 RTP PGA
hold on
contourf(x, y, map_2475, [0:0.2:3.2])
plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
hold off
caxis([0, 3.0])
cb = colorbar("FontName", font, "Fontsize", tick_fs);
cb.Label.String = "a_{g} [m/s²]";
colormap(p2, "jet")
ax = gca;
ax.XAxis.FontSize = tick_fs;
ax.XAxis.FontName = font;
ax.YAxis.FontSize = tick_fs;
ax.YAxis.FontName = font;
set(ax,'Units','centimeters','Position', pos2, "YTickLabel", [])
set(get(p2, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)

annotation(textbox=[0.6, 0.58, 0.3, 0.3], ...
           String=("Max. a_{g,2475} = " + string(round(PGA_2475_median_max,3)) + " m·s^{-2}" + newline + ...
                   "Min. a_{g,2475} = " + string(round(PGA_2475_median_min,3)) + " m·s^{-2}"), ...
           FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs, Color="white");

annotation(textbox=[0.515, 0.88, 0.1, 0.1], String="(b)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

saveas(f, figure_folder + "\" + save_file + ".png")
saveas(f, figure_folder + "\" + save_file + ".pdf")
