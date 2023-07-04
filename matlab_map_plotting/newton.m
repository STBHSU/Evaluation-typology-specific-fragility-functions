function [x2] = newton(x1, fx1, x0, fx0)
%NEWTON applies newtons method to produce the next x value
%   This function uses the slope between the current point (x1, y1) and the
%   previous point (x0, y0) as an approximation for the tanget at (x1, y1).

fx1dx = (fx1 - fx0) / (x1 - x0);
x2 = x1 - (fx1 / fx1dx);

% fprintf("fx1dx: %f\n", fx1dx)

end