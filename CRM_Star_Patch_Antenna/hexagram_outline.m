function pts = hexagram_outline(Rtip)
% hexagram_outline  Centre-line of a regular 6-point star (Star of David).
%
%   pts = hexagram_outline(Rtip)
%
% Returns a closed 2x13 polyline: 12 vertices alternating between tip
% radius Rtip and inner radius Rtip/sqrt(3), first tip pointing up (+y),
% matching the paper figures.
Ri = Rtip / sqrt(3);
n = 12;
pts = zeros(2, n + 1);
for k = 1:n
    ang = (90 - (k-1)*30) * pi/180;   % tip at 90 deg, then every 30 deg
    if mod(k-1, 2) == 0
        r = Rtip;                     % tip vertex
    else
        r = Ri;                       % inner (notch) vertex
    end
    pts(:, k) = [r*cos(ang); r*sin(ang)];
end
pts(:, n+1) = pts(:, 1);              % close the loop
end
