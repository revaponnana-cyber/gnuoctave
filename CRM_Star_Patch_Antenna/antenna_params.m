function p = antenna_params()
% antenna_params  Design parameters for the CRM-inspired star-shaped
% microstrip patch antenna (FR-4, 50x50x1.6 mm^3, f0 = 2.82 GHz).
%
% The values below follow Table 1 (antenna) and Table 2 (CRM unit cell)
% of the paper. Figure-derived / assumed values are marked ASSUMPTION and
% can be tuned in seconds; run preview_layout to see the effect.
%
% Units: all lengths in mm, frequencies in Hz.
%
% Paper mapping (Table 1):
%   A     = 14.5  patch radius (confirmed by paper text, Sect. 4.1.2)
%   B,D,E = 4.6 / 4.04 / 1.84  star "tip-to-notch" horizontal lengths
%           (Fig. 4a) -> tip radii 2x = [9.2, 8.08, 3.68]
%   C     = 0.75  outer ring-slot width (Fig. 4a)
%   F     = 1.24  star outline slot width (Fig. 4a, inner star detail)
%   G     = 2.43  ring-slot opening width at bottom (ASSUMPTION)
%   H,I   = 6 / 6.66 bottom-point extents; used as cross-check only
%   J     = 20    DGS cross arm TOTAL length (Fig. 4b: arms span 20 mm)
%   K     = 1.2   DGS slot width (Fig. 4b)
%   L     = 1     coax pin RADIUS -> 2 mm diameter dot (Fig. 4b scale)
% Paper mapping (Table 2, CRM unit cell, Fig. 5):
%   T,U   = 5 / 5.95 inner/outer ring radii (centerlines)
%   M     = 0.2   ring trace width (-> 0.25 after mesh-grid quantize)
%   Q     = 0.3   ring split-gap width (-> 0.25 after quantize)
%   R     = 1.5   inner disk DIAMETER (Fig. 5 scale)
%   N     = 0.8   outer corner disk DIAMETER (Fig. 5 scale)
%   S     = 0.95  stub length, disk edge to inner-ring edge
%   P     = 0.44  stub trace width (-> 0.5 after quantize)
%   O     = 0.95  ring pitch U-T (descriptive, not a construction dim)

p.unit = 1e-3;                 % mm -> m

% ---------- substrate / design frequency ----------
p.sub.L   = 50;                % Sub-L (x)
p.sub.W   = 50;                % Sub-W (y)
p.sub.h   = 1.6;               % thickness
p.sub.epsR = 4.4;              % FR-4
p.sub.tanD = 0.025;            % FR-4 loss tangent
p.f_design = 2.82e9;           % design / resonance target

% ---------- circular patch + slots (top metal, z = h) ----------
p.patch.A = 14.5;                                % patch radius
p.star.tipR  = 2 * [4.6, 4.04, 1.84];            % outer/mid/inner tip radii
p.star.slotW = [0.75, 0.75, 1.24];               % outline widths [C, C, F];
                                                 % outer/mid use C so the
                                                 % nested slots stay separate
p.ring.slotW  = 0.75;          % C
p.ring.rOuter = 13.5;          % ASSUMPTION: ring outer edge 1 mm inside patch
p.ring.gapW   = 2.43;          % G, opening at bottom (ASSUMPTION)

% ---------- CRM unit cells (top metal, 4x corners) ----------
p.crm.enabled = true;          % set false to reproduce "without CRM" case
p.crm.pos     = 16.8;          % |x|,|y| of cell centres (ASSUMPTION, Fig.2a)
p.crm.rIn     = 5;             % T, inner-ring centreline radius
p.crm.rOut    = 5.95;          % U, outer-ring centreline radius
p.crm.trace   = 0.2;           % M (quantized to mesh grid in model)
p.crm.gap     = 0.3;           % Q (quantized to mesh grid in model)
p.crm.diskD   = 1.5;           % R, inner disk diameter
p.crm.stubLen = 0.95;          % S, disk edge -> inner-ring edge
p.crm.stubW   = 0.44;          % P (quantized to mesh grid in model)
p.crm.outDiskD = 0.8;          % N, outer disk diameter
p.crm.outDiskR = 7.6;          % outer-disk centre radius (ASSUMPTION, Fig.5)

% ---------- cross-shaped DGS (ground plane, z = 0) ----------
p.dgs.enabled = true;          % set false for a solid ground plane
p.dgs.arm     = 20;            % J, total arm length (slot spans +/-10)
p.dgs.w       = 1.2;           % K, slot width
p.dgs.center  = [0, 14];       % cross centre (ASSUMPTION, Fig. 4b scale)

% ---------- coaxial probe feed (modelled as 50 Ohm lumped port) ----------
% Fig. 4b places the feed dot at about (0,-3.5), but with the Table 1 star
% sizes that point lands on the inner-star slot. It is shifted down the
% centre line to (0,-5.5) so the probe contacts patch metal; tune for the
% 50 Ohm match (typical circle-patch match radius is 4-6 mm here).
p.feed.pos = [0, -5.5];        % (x,y), figure-derived + slot clearance
p.feed.R   = 50;               % coax / port impedance
p.feed.pinR = 1;               % L, inner-pin radius (drawing only)

% ---------- FDTD excitation / run control ----------
p.f0 = 2.5e9;                  % Gaussian pulse centre (covers 1-4 GHz)
p.fc = 1.5e9;                  % 20 dB corner (fmax = f0+fc = 4 GHz)
p.NrTs = 50000;                % max time steps
p.EndCriteria = 1e-4;          % stop when energy decayed by 1e-4
p.SimBox = [110, 110, 90];     % full simulation box (MUR boundaries)
p.mesh_snap = 0.25;            % quantize staircase geometry to this grid (mm)

% ---------- outputs ----------
p.Sim_Path = 'tmp_CRM_Star_Patch';
p.Sim_CSX  = 'crm_star_patch.xml';
p.run_simulation = true;       % false -> write XML + plot geometry only
p.show_geometry  = true;       % CSXGeomPlot preview before running
end
