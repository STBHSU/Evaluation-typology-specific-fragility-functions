function [map_mat, x, y] = matrix_for_contours(data, lat, lon, dlon, dlat)
%MATRIX_FOR_CONTOURS distributes a vector of points into a matrix 
%   Used to create contour maps from the vector of data that is related by
%   latitude and longitude vectors with constant steps defined by dlon and
%   dlat

% create 2D matrix of points for hazard contour plot
maxcoord = length(lat);
y = min(lat):dlat:max(lat);
x = min(lon):dlon:max(lon);
map_mat = zeros(max(size(y)),max(size(x)));

for coord=1:maxcoord
        p = 1+ round((lon(coord)-min(lon))/dlon);
        q = 1+ round((lat(coord)-min(lat))/dlat);
        
        map_mat(q,p) = data(coord);
end

map_mat(map_mat == 0) = NaN;
end