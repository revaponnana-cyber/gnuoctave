%% preview_layout  2D top/bottom layout of the CRM star patch antenna.
%% Pure MATLAB/Octave (no openEMS calls): run this first to verify every
%% dimension before launching the full-wave simulation.
close all
p = antenna_params();
h = p.sub.h;

figure('Name', 'CRM star patch layout', 'Position', [80, 80, 1150, 520]);

%% ---------- top view ----------
subplot(1, 2, 1);
hold on; axis equal; grid on;
xlim([-26, 26]); ylim([-26, 26]);
title('Top: star-slotted patch + 4x CRM cells');
xlabel('x / mm'); ylabel('y / mm');
rectangle('Position', [-25, -25, 50, 50], 'FaceColor', [0.93, 0.93, 0.93], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);                       % substrate

th = linspace(0, 2*pi, 361);                                   % patch disk
patch(p.patch.A*cos(th), p.patch.A*sin(th), [0.85, 0.15, 0.15], ...
    'EdgeColor', 'k');

for s = 1:3                                                   % star slots
    pts = hexagram_outline(p.star.tipR(s));
    plot(pts(1, :), pts(2, :), 'w-', 'LineWidth', 3);
end
rOut = p.ring.rOuter; rIn = rOut - p.ring.slotW;              % ring slot
gapHalf = (p.ring.gapW/2) / (rOut - p.ring.slotW/2);
a = linspace(-pi/2 + gapHalf, 3*pi/2 - gapHalf, 200);
plot(rOut*cos(a), rOut*sin(a), 'w-', 'LineWidth', 2);
plot(rIn*cos(a), rIn*sin(a), 'w-', 'LineWidth', 2);

for sx = [-1, 1]                                              % CRM cells
    for sy = [-1, 1]
        draw_crm_cell(sx*p.crm.pos, sy*p.crm.pos, p, th);
    end
end
plot(p.feed.pos(1), p.feed.pos(2), 'b+', 'MarkerSize', 10, 'LineWidth', 2);

%% ---------- bottom view ----------
subplot(1, 2, 2);
hold on; axis equal; grid on;
xlim([-26, 26]); ylim([-26, 26]);
title('Bottom: cross DGS + probe feed');
xlabel('x / mm'); ylabel('y / mm');
rectangle('Position', [-25, -25, 50, 50], 'FaceColor', [0.72, 0.45, 0.35], ...
    'EdgeColor', 'k', 'LineWidth', 1.5);                       % ground
cx0 = p.dgs.center(1); cy0 = p.dgs.center(2);
ax = p.dgs.arm/2; w = p.dgs.w/2;                              % cross slot
rectangle('Position', [cx0-ax, cy0-w, 2*ax, 2*w], 'FaceColor', 'w', ...
    'EdgeColor', 'k');
rectangle('Position', [cx0-w, cy0-ax, 2*w, 2*ax], 'FaceColor', 'w', ...
    'EdgeColor', 'k');
rectangle('Position', ...                                     % feed pin
    [p.feed.pos(1)-p.feed.pinR, p.feed.pos(2)-p.feed.pinR, ...
     2*p.feed.pinR, 2*p.feed.pinR], ...
    'Curvature', [1, 1], 'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);
drawnow;

function draw_crm_cell(cx, cy, p, th)
% Scaled drawing of one CRM unit cell (Fig. 5), metal in dark red.
rIn = p.crm.rIn; rOut = p.crm.rOut; tr = max(p.crm.trace, 0.2);
col = [0.75, 0.1, 0.1];
plot(cx + rIn*cos(th), cy + rIn*sin(th), 'Color', col, 'LineWidth', 2.5);
plot(cx + rOut*cos(th), cy + rOut*sin(th), 'Color', col, 'LineWidth', 2.5);
diskR = p.crm.diskD/2;
diskRp = p.crm.rIn - tr/2 - p.crm.stubLen - diskR;
for ang = [0, 90, 180, 270]
    dx = cosd(ang); dy = sind(ang);
    px = cx + diskRp*dx; py = cy + diskRp*dy;
    patch(px + diskR*cos(th), py + diskR*sin(th), col, 'EdgeColor', col);
    if abs(dx) > 0.5
        rectangle('Position', [min(cx, px), py - p.crm.stubW/2, ...
            abs(px - cx), p.crm.stubW], 'FaceColor', col, 'EdgeColor', col);
    else
        rectangle('Position', [px - p.crm.stubW/2, min(cy, py), ...
            p.crm.stubW, abs(py - cy)], 'FaceColor', col, 'EdgeColor', col);
    end
end
oR = p.crm.outDiskD/2; oOff = p.crm.outDiskR/sqrt(2);
for qx = [-1, 1]
    for qy = [-1, 1]
        px = cx + qx*oOff; py = cy + qy*oOff;
        patch(px + oR*cos(th), py + oR*sin(th), col, 'EdgeColor', col);
    end
end
end
