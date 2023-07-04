function [theta] = theta_linear(m, c, im_des)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
theta = m .* im_des + c;
end