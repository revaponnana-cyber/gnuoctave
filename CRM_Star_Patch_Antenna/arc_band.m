function pts = arc_band(rOuter, rInner, a0deg, a1deg, n)
% arc_band  Filled annular-sector polygon (for CRM split-ring arcs).
%
%   pts = arc_band(rOuter, rInner, a0deg, a1deg, n)
%
% Angles in degrees, swept counter-clockwise with n points per arc.
% Returns a closed 2x(2n+2) polygon: outer arc out, inner arc back.
if nargin < 5
    n = 16;
end
a = linspace(a0deg, a1deg, n) * pi/180;
outer = [rOuter*cos(a); rOuter*sin(a)];
inner = [rInner*cos(fliplr(a)); rInner*sin(fliplr(a))];
pts = [outer, inner];
pts = [pts, pts(:, 1)];               % close the loop
end
