function vq = snap_coords(v, s)
% snap_coords  Quantize coordinates to a uniform grid of step s.
% Keeps the FDTD mesh-line count bounded when many staircase boxes share
% edges. Works element-wise on any size array.
vq = round(v / s) * s;
end
