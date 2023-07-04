function  [] = plot_rt_pga_maps(save_name, fig_folder, lat, lon, Cr, pga_risk, dlat, dlon, figure_name)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% INPUTS

% file containing variables for mapping
% map_file = "map_values_steel_mrfs.mat";
% save_name = "Fig16_pga-risk-map_steel-mrfs";

% fig_folder = "C:\Users\Nicholas Clemett\Documents\risiko-basierte-erdbebenkarte\Bearbeitung\figures";

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

%% PLOTTING

% load(map_file)

% create 2D matrix of points for contour plot
maxcoord = length(lat);
% dlat =0.1;
% dlon =0.1;
y = min(lat):dlat:max(lat);
x = min(lon):dlon:max(lon);
Cr_map = zeros(max(size(y)),max(size(x)));
pga_risk_map = zeros(max(size(y)),max(size(x)));

Cr_max = max(Cr);
Cr_min = min(Cr);
PGA_max = max(pga_risk);
PGA_min = min(pga_risk);

for coord=1:maxcoord
        p = 1+ round((lon(coord)-min(lon))/dlon);
        q = 1+ round((lat(coord)-min(lat))/dlat);
        
        Cr_map(q,p) = Cr(coord);
        pga_risk_map (q,p) = pga_risk(coord);
end

Cr_map(Cr_map == 0) = NaN;
pga_risk_map(pga_risk_map==0)=NaN;

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
f=figure("Units","centimeters", "Position",[10,10,fig_w,fig_h], ...
         "Name", figure_name, "NumberTitle", "off");

p1 = subplot(1,2,1); % first subplot for risk-based pga
hold on
contourf(x, y, pga_risk_map,[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0])%, edgecolor='none')
plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
hold off
cb = colorbar("FontName", font, "Fontsize", tick_fs); 
cb.Label.String = "a_{g,risk} [m/s²]";
colormap(p1, "jet")
caxis([0, 2.0])
ax = gca;
set(ax,'Units','centimeters', "Position", pos1)
ax.XAxis.FontSize = tick_fs;
ax.XAxis.FontName = font;
ax.YAxis.FontSize = tick_fs;
ax.YAxis.FontName = font;
set(get(p1, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)
set(get(p1, "Ylabel"), "String", "Latitude (°)", "Fontsize", label_fs, "Fontname", font)

annotation(textbox=[0.15, 0.58, 0.3, 0.3], ...
           String=("Max. a_{g,risk} = " + num2str(PGA_max, '%.3f') + " m·s^{-2}" + newline + ...
                   "Min. a_{g,risk} = " + num2str(PGA_min, '%.3f') + " m·s^{-2}"), ...
           FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);

annotation(textbox=[0.09, 0.88, 0.1, 0.1], String="(a)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

mymap = redgreyblue_cmap2(1, 3, 5, 20);

p2 = subplot(1,2,2); % second subplot is the risk-coefficient
hold on
contourf(x, y, Cr_map, [0:0.2:14], edgecolor='none')
plot(border(:,1), border(:,2), "LineWidth",1.5, "Color","black")
hold off
caxis([0.8, 5])
cb = colorbar("FontName", font, "Fontsize", tick_fs);
cb.Label.String = "Risk Coefficient";
colormap(p2,mymap)
% colormap(p2, "copper")
% colormap(p2, flipud(colormap(p2)))
ax = gca;
ax.XAxis.FontSize = tick_fs;
ax.XAxis.FontName = font;
ax.YAxis.FontSize = tick_fs;
ax.YAxis.FontName = font;
set(ax,'Units','centimeters','Position', pos2, "YTickLabel", [])
set(get(p2, "Xlabel"), "String", "Longitude (°)", "Fontsize", label_fs, "Fontname", font)

annotation(textbox=[0.65, 0.55, 0.3, 0.3], ...
           String=("Max. CR = " + num2str(Cr_max, '%.3f') + newline + ...
                   "Min. CR = " + num2str(Cr_min, '%.3f')), ...
           FitBoxToText="on", EdgeColor="none", FontName=font, FontSize=tick_fs);

annotation(textbox=[0.55, 0.88, 0.1, 0.1], String="(b)", EdgeColor="none", ...
        FitBoxToText="on", FontName=font, FontSize=tick_fs)

saveas(f, fig_folder + "\" + save_name + ".png")
saveas(f, fig_folder + "\" + save_name + ".pdf")

end