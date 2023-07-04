% Title: Probability of Collapse Maps 
% Author: Nicholas Clemett
% Date: 14.09.22

% Description:
%   - plots the maps of germany showing the spatial distribution of the
%     probability of collapse --> Figure 14 in paper

clear 
close all
clc

%% Input Parameters
% load and save file names
key = "fig14_pc_";
key_2018 = "_NA-2018";
figure_label = "Pc";
data_folder = "data_out";
figure_folder = pwd;

% list of indentifiers of fragility curves to use
% ids = ["rc-mrf-m-pw"];
ids = ["rc-mrf-m-rto"];
% ids = ["rc-mrf-r" "rc-mrf-m" "s-mrf" "rc-wds"];
% ids = ["rc-mrf-r" "rc-mrf-m" "s-mrf"];

% constants
dlat = 0.1;
dlon = 0.1;

% figure size, aspect ratio im_cstand offsets
fig_w = 16;
img_ar = 0.83; 
lft = 1.5; 
gap = 2.2;
rgt = 2.0;  
bot = 1.0;  
top = 0.2;  

% Font Settings
font = "Times";
label_fs = 12;
tick_fs = 9;

% Create a new figure for the maps
% image and figure sizes
img_w = (fig_w - (lft + gap + rgt)) / 2;
img_h = img_w / img_ar;
fig_h = img_h + bot + top;

% position vectors for two images side by side
pos1 = [lft, bot, img_w, img_h];
pos2 = [lft + img_w + gap, bot, img_w, img_h];

% Create a new figure for the maps
figure_name = "FigX_pc-RC-MRF-M";

%% Plotting
% load data and plot
for ii = 1:1:length(ids)
    data = load(fullfile(pwd, data_folder, key + ids(ii) + key_2018 + ".mat"));
    [map, x, y] = matrix_for_contours(data.pc, data.lat, data.lon, dlat, dlon);

    map_max = max(max(map));
    map_min = min(min(map));


    
    % import border points
    border = readmatrix("new_border.csv");
    
    % Create Figure
    f=figure("Units","centimeters", "Position",[10,10,fig_w,fig_h]);
    p1 = subplot(1,2,1); % first subplot for risk-based pga
    contourf(x, y, map, [0:0.1:1.8])
    hold on
    plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
    hold off
    caxis([0, 1.8])
    cb = colorbar("FontName", font, "Fontsize", tick_fs);
    cb.Label.String = "P[C] in 50 years (%)";
    colormap("jet")
    ax = gca;
    set(ax,'Units','centimeters', "Position", pos1)
    ax.XAxis.FontSize = tick_fs;
    ax.XAxis.FontName = font;
    ax.YAxis.FontSize = tick_fs;
    ax.YAxis.FontName = font;
    set(get(p1, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
    set(get(p1, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

%     title(figure_label + "-" + ids(ii))
    
    % adding max and minimum values
    annotation(textbox=[0.18, 0.55, 0.3, 0.3], ...
               String=("Max. = " + string(round(map_max,3)) + "%" + newline + ...
                       "Min. = " + string(round(map_min,3)) + "%"), ...
               FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);
    
    annotation(textbox=[0.09, 0.88, 0.1, 0.1], String="(a)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

%     saveas(f, fullfile(figure_folder, figure_name + ".png"))
%     saveas(f, fullfile(figure_folder, figure_name + ".pdf"))

    
end

ids = ["rc-mrf-m-pw"];

for ii = 1:1:length(ids)
    data = load(fullfile(pwd, data_folder, key + ids(ii) + key_2018 + ".mat"));
    [map, x, y] = matrix_for_contours(data.pc, data.lat, data.lon, dlat, dlon);

    map_max = max(max(map));
    map_min = min(min(map));
    
    % import border points
    border = readmatrix("new_border.csv");
    
    % Create Figure
    p2 = subplot(1,2,2); % first subplot for risk-based pga
    contourf(x, y, map, [0:0.1:1.8])
    hold on
    plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
    hold off
    caxis([0, 1.8])
    cb = colorbar("FontName", font, "Fontsize", tick_fs);
    cb.Label.String = "P[C] in 50 years (%)";
    colormap("jet")
    ax = gca;
    set(ax,'Units','centimeters', "Position", pos2, "YTickLabel", [])
    ax.XAxis.FontSize = tick_fs;
    ax.XAxis.FontName = font;
    ax.YAxis.FontSize = tick_fs;
    ax.YAxis.FontName = font;
    set(get(p2, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
%     set(get(p2, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

%     title(figure_label + "-" + ids(ii))
    
    % adding max and minimum values
    annotation(textbox=[0.65, 0.55, 0.3, 0.3], ...
               String=("Max. = " + string(round(map_max,3)) + "%" + newline + ...
                       "Min. = " + sprintf('%.3f', round(map_min,3)) + "%"), ...
               FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);
    
    annotation(textbox=[0.55, 0.88, 0.1, 0.1], String="(b)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

    saveas(f, fullfile(figure_folder, figure_name + ".png"))
    saveas(f, fullfile(figure_folder, figure_name + ".pdf"))

    
end




