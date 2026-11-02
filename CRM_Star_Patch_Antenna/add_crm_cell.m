function CSX = add_crm_cell(CSX, p, cx, cy, zTop)
% add_crm_cell  Place one circular-ring-metamaterial unit cell (Fig. 5).
%
% Two concentric split rings (gaps at cardinals/diagonals), four inner
% disks joined by stubs, four outer corner disks joined by diagonal stubs.
% Ring radii/traces are quantized to the mesh grid (<=0.13 mm shift).
snap = p.mesh_snap;

% Ring bands on the mesh grid (trace 0.2 -> 0.25, pitch stays ~1.0 = O).
rInA  = snap_coords(p.crm.rIn  - p.crm.trace/2, snap);   % 4.75
rInB  = snap_coords(p.crm.rIn  + p.crm.trace/2, snap);   % 5.00
rOutA = snap_coords(p.crm.rOut - p.crm.trace/2, snap);   % 5.75
rOutB = snap_coords(p.crm.rOut + p.crm.trace/2, snap);   % 6.00
gapHalfIn  = (snap_coords(p.crm.gap, snap)/2) / p.crm.rIn  * 180/pi;
gapHalfOut = (snap_coords(p.crm.gap, snap)/2) / p.crm.rOut * 180/pi;

% Inner ring: 4 arcs with split gaps at 0/90/180/270 deg.
for q = 0:3
    a0 = q*90 + 45 + gapHalfIn;
    a1 = q*90 + 135 - gapHalfIn;
    pts = arc_band(rInB, rInA, a0, a1, 14);
    pts = snap_polygon(pts + repmat([cx; cy], 1, size(pts, 2)), snap);
    CSX = AddPolygon(CSX, 'crm', 10, 'z', zTop, pts);
end

% Outer ring: 4 arcs with split gaps at 45/135/225/315 deg.
for q = 0:3
    a0 = q*90 + gapHalfOut;
    a1 = q*90 + 90 - gapHalfOut;
    pts = arc_band(rOutB, rOutA, a0, a1, 14);
    pts = snap_polygon(pts + repmat([cx; cy], 1, size(pts, 2)), snap);
    CSX = AddPolygon(CSX, 'crm', 10, 'z', zTop, pts);
end

% Inner disks (diameter R) + cardinal stubs (width P -> grid).
diskR  = p.crm.diskD/2;
diskRp = snap_coords(p.crm.rIn - p.crm.trace/2 - p.crm.stubLen - diskR, snap);
stubW  = max(snap_coords(p.crm.stubW, snap), snap);
stubIn = diskRp - diskR + 0.25;             % overlap into the disk
stubOut = rOutA - stubW;                    % stop ~one width before outer ring
for ang = [0, 90, 180, 270]
    dx = cos(ang*pi/180); dy = sin(ang*pi/180);
    pts = circle_points(diskR, 20, cx + diskRp*dx, cy + diskRp*dy);
    pts = snap_polygon(pts, snap);
    CSX = AddPolygon(CSX, 'crm', 10, 'z', zTop, pts);
    if abs(dx) > 0.5                            % horizontal stub (box)
        xa = snap_coords(min(cx + stubIn*dx, cx + stubOut*dx), snap);
        xb = snap_coords(max(cx + stubIn*dx, cx + stubOut*dx), snap);
        ya = snap_coords(cy - stubW/2, snap);
        yb = snap_coords(cy + stubW/2, snap);
        CSX = AddBox(CSX, 'crm', 10, [xa, ya, zTop], [xb, yb, zTop]);
    else                                        % vertical stub (box)
        xa = snap_coords(cx - stubW/2, snap);
        xb = snap_coords(cx + stubW/2, snap);
        ya = snap_coords(min(cy + stubIn*dy, cy + stubOut*dy), snap);
        yb = snap_coords(max(cy + stubIn*dy, cy + stubOut*dy), snap);
        CSX = AddBox(CSX, 'crm', 10, [xa, ya, zTop], [xb, yb, zTop]);
    end
end

% Outer corner disks (diameter N) + diagonal stub staircases.
oDiskR = p.crm.outDiskD/2;
oOff   = snap_coords(p.crm.outDiskR/sqrt(2), snap);
for sx = [-1, 1]
    for sy = [-1, 1]
        ox = cx + sx*oOff; oy = cy + sy*oOff;
        pts = circle_points(oDiskR, 12, ox, oy);
        pts = snap_polygon(pts, snap);
        CSX = AddPolygon(CSX, 'crm', 10, 'z', zTop, pts);
        ux = sx/sqrt(2); uy = sy/sqrt(2);       % diagonal unit vector
        rA = p.crm.rOut + p.crm.trace/2 - 0.25; % overlap outer ring
        rB = norm([ox, oy] - [cx, cy]) - oDiskR + 0.25;
        seg = [cx + rA*ux, cx + rB*ux; cy + rA*uy, cy + rB*uy];
        CSX = staircase_boxes(CSX, 'crm', 10, seg, stubW, zTop, zTop, snap);
    end
end
end
