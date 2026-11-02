function CSX = staircase_boxes(CSX, prop, prio, pts, w, z0, z1, snap)
% staircase_boxes  Approximate a polyline band with small axis-aligned boxes.
%
%   CSX = staircase_boxes(CSX, prop, prio, pts, w, z0, z1, snap)
%
% Used for diagonal slot edges (star outlines, ring slot -> air boxes that
% cut the patch metal by higher priority) and diagonal metal stubs.
%   pts  : 2xN polyline (star centrelines, circle paths, stub paths)
%   w    : band width; boxes are w x w, stepped at w/2 for full overlap
%   z0,z1: box z-extent (use z0 == z1 for a zero-thickness metal sheet)
%   snap : quantize every box edge to this grid to bound the mesh count
ws = max(snap_coords(w, snap), snap);
step = max(ws/2, snap);
for s = 1:(size(pts, 2)-1)
    p0 = pts(:, s);
    p1 = pts(:, s+1);
    segLen = norm(p1 - p0);
    if segLen == 0
        continue;
    end
    nStep = max(1, ceil(segLen / step));
    for k = 0:nStep
        c = p0 + (p1 - p0) * (k / nStep);
        x0 = snap_coords(c(1) - ws/2, snap);
        x1 = snap_coords(c(1) + ws/2, snap);
        y0 = snap_coords(c(2) - ws/2, snap);
        y1 = snap_coords(c(2) + ws/2, snap);
        if x1 <= x0
            x1 = x0 + snap;
        end
        if y1 <= y0
            y1 = y0 + snap;
        end
        CSX = AddBox(CSX, prop, prio, [x0, y0, z0], [x1, y1, z1]);
    end
end
end
