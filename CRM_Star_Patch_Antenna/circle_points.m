function pts = circle_points(r, n, cx, cy)
% circle_points  Closed 2x(n+1) polygon approximating a circle of radius r
% centred at (cx, cy). Used for the patch disk and CRM disks.
if nargin < 2 || isempty(n)
    n = 72;
end
if nargin < 4
    cx = 0; cy = 0;
end
a = linspace(0, 2*pi, n+1);
pts = [cx + r*cos(a); cy + r*sin(a)];
end
