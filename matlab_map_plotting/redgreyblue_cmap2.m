function [cmap] = redgreyblue_cmap2(min, mid, max, steps)
% cmap = redgreyblue_cmap2(0, 1, 2, 10)
% redgreyblue_cmap colormap from red through grey to blue
%   produces a colour map with the specified number of steps either side of
%   the middle grey colour.

% interpolate between red = [1 0 0]; blue = [0 0 1]; grey = [200, 200, 200]
if ((max - mid) < (mid -min))
    u_steps = steps;
    d_steps = steps * (mid -min) / (max - mid);

elseif ((max - mid) == (mid -min))
    u_steps = steps;
    d_steps = steps;
else
    d_steps = steps;
    u_steps = steps * (max - mid) / (mid -min);
end

main_channel_u = linspace(1, 230/255, u_steps)';
other_channel_u = linspace(0, 230/255, u_steps)';

main_channel_d = linspace(1, 230/255, d_steps)';
other_channel_d = linspace(0, 230/255, d_steps)';

cmap = [other_channel_d other_channel_d main_channel_d; 
        flip(main_channel_u(1:(end-1))) flip(other_channel_u(1:(end-1))) flip(other_channel_u(1:(end-1)))];

end