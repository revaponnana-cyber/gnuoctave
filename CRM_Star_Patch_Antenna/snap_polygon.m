function pts = snap_polygon(pts, s)
% snap_polygon  Snap a 2xN polygon/polyline to grid s and drop consecutive
% duplicate points (which CSXCAD would reject).
pts = snap_coords(pts, s);
keep = true(1, size(pts, 2));
for k = 2:size(pts, 2)
    if all(pts(:, k) == pts(:, k-1))
        keep(k) = false;
    end
end
pts = pts(:, keep);
end
