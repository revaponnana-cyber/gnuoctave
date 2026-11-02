%% CRM-inspired star-shaped microstrip patch antenna (2.82 GHz, FR-4)
%% ------------------------------------------------------------------
%% Full-wave openEMS model of the paper design: circular patch with triple
%% star slots + outer ring slot, 4x corner CRM unit cells, cross-shaped DGS
%% and a 50 Ohm coaxial-probe feed (modelled as a lumped port).
%%
%% Runs in MATLAB (R2016b or newer) or GNU Octave with the openEMS +
%% CSXCAD matlab/octave interface on the path. Geometry helpers are plain
%% m-code, so preview_layout.m and circular_patch_theory.m run anywhere.
%%
%% Paper targets (simulated): S11 = -25.81 dB @ 2.83 GHz, gain 7.3 dB,
%% directivity 6.6 dB, radiation efficiency 96 %.
%% NOTE: gain cannot physically exceed directivity (G = eta*D), so treat
%% the two reported peak numbers with care; this script computes both from
%% the same NF2FF solution. See README.md.

close all
clear
clc

%% 0. Parameters (paper Tables 1-2, see antenna_params for the mapping)
p = antenna_params();
physical_constants;             % EPS0, MU0, C0, c0 (openEMS)

unit = p.unit;                  % all lengths in mm
h  = p.sub.h;                   % substrate thickness
Lx = p.sub.L; Ly = p.sub.W;

% Dielectric loss from the FR-4 loss tangent at the design frequency.
kappa = 2*pi*p.f_design*EPS0*p.sub.epsR*p.sub.tanD;

%% 1. FDTD settings: broadband pulse covering the 1-4 GHz band
FDTD = InitFDTD('NrTS', p.NrTs, 'EndCriteria', p.EndCriteria);
FDTD = SetGaussExcite(FDTD, p.f0, p.fc);
FDTD = SetBoundaryCond(FDTD, {'MUR','MUR','MUR','MUR','MUR','MUR'});

%% 2. Materials: substrate + "air" used to cut slots by higher priority
CSX = InitCSX();
CSX = AddMaterial(CSX, 'substrate');
CSX = SetMaterialProperty(CSX, 'substrate', 'Epsilon', p.sub.epsR, ...
                                                 'Kappa', kappa);
CSX = AddMaterial(CSX, 'slot_air');
CSX = SetMaterialProperty(CSX, 'slot_air', 'Epsilon', 1, 'Kappa', 0);
CSX = AddMetal(CSX, 'patch');
CSX = AddMetal(CSX, 'crm');
CSX = AddMetal(CSX, 'gnd');

CSX = AddBox(CSX, 'substrate', 0, [-Lx/2, -Ly/2, 0], [Lx/2, Ly/2, h]);

%% 3. Circular patch + triple star slots + outer ring slot (z = h)
snap = p.mesh_snap;
patchPts = snap_polygon(circle_points(p.patch.A, 72), snap);
CSX = AddPolygon(CSX, 'patch', 10, 'z', h, patchPts);

% Slot air boxes pierce the patch sheet (priority 20 > 10).
zS0 = h - 0.2; zS1 = h + 0.2;
for s = 1:3
    pts = snap_polygon(hexagram_outline(p.star.tipR(s)), snap);
    CSX = staircase_boxes(CSX, 'slot_air', 20, pts, ...
                          p.star.slotW(s), zS0, zS1, snap);
end

% Outer ring slot (annular path, opening of width G left at the bottom).
rC = p.ring.rOuter - p.ring.slotW/2;            % centreline radius
gapHalf = (p.ring.gapW/2) / rC;                 % half-angle of opening
nRing = 160;
aRing = linspace(-pi/2 + gapHalf, 3*pi/2 - gapHalf, nRing);
ringPts = snap_polygon([rC*cos(aRing); rC*sin(aRing)], snap);
CSX = staircase_boxes(CSX, 'slot_air', 20, ringPts, ...
                      p.ring.slotW, zS0, zS1, snap);

%% 4. CRM unit cells at the four substrate corners (z = h)
if p.crm.enabled
    for sx = [-1, 1]
        for sy = [-1, 1]
            CSX = add_crm_cell(CSX, p, sx*p.crm.pos, sy*p.crm.pos, h);
        end
    end
else
    disp('CRM cells disabled (ablation case "without CRM").');
end

%% 5. Ground plane + cross-shaped DGS (z = 0)
CSX = AddBox(CSX, 'gnd', 10, [-Lx/2, -Ly/2, 0], [Lx/2, Ly/2, 0]);
if p.dgs.enabled
    cx0 = p.dgs.center(1); cy0 = p.dgs.center(2);
    ax = p.dgs.arm/2; w = snap_coords(p.dgs.w, snap)/2;
    zG0 = -0.5; zG1 = 0.2;
    CSX = AddBox(CSX, 'slot_air', 20, ...
        [snap_coords(cx0-ax, snap), snap_coords(cy0-w, snap), zG0], ...
        [snap_coords(cx0+ax, snap), snap_coords(cy0+w, snap), zG1]);
    CSX = AddBox(CSX, 'slot_air', 20, ...
        [snap_coords(cx0-w, snap), snap_coords(cy0-ax, snap), zG0], ...
        [snap_coords(cx0+w, snap), snap_coords(cy0+ax, snap), zG1]);
else
    disp('DGS disabled (solid ground plane).');
end

%% 6. Coaxial probe feed as a 50 Ohm vertical lumped port
fx = p.feed.pos(1); fy = p.feed.pos(2);
[CSX, port] = AddLumpedPort(CSX, 5, 1, p.feed.R, ...
    [fx, fy, 0], [fx, fy, h], [0, 0, 1], true);

%% 7. Mesh: detect all edges, refine metal edges, smooth to lambda/20
mesh.x = [-p.SimBox(1)/2, p.SimBox(1)/2];
mesh.y = [-p.SimBox(2)/2, p.SimBox(2)/2];
mesh.z = unique([-p.SimBox(3)/3, 2*p.SimBox(3)/3, -12, -4, -0.5, 0, ...
                 linspace(0, h, 7), h+0.2, h+0.5, 3, 5, 8, 12, 18, ...
                 26, 36, 48]);
mesh = DetectEdges(CSX, mesh);
mesh = DetectEdges(CSX, mesh, 'SetProperty', 'patch', ...
    '2D_Metal_Edge_Res', c0/(p.f0+p.fc)/unit/50);
mesh = DetectEdges(CSX, mesh, 'SetProperty', 'crm', ...
    '2D_Metal_Edge_Res', c0/(p.f0+p.fc)/unit/50);
mesh = SmoothMesh(mesh, c0/(p.f0+p.fc)/unit/20);
mesh.x = snap_coords(mesh.x, snap/5);
mesh.y = snap_coords(mesh.y, snap/5);
mesh.x = unique(mesh.x); mesh.y = unique(mesh.y); mesh.z = unique(mesh.z);
CSX = DefineRectGrid(CSX, unit, mesh);
disp(['Mesh lines: nx=' num2str(length(mesh.x)) ...
      ', ny=' num2str(length(mesh.y)) ...
      ', nz=' num2str(length(mesh.z)) ...
      ' -> cells~' num2str(length(mesh.x)*length(mesh.y)*length(mesh.z))]);

% H-field dump at resonance for surface-current inspection (paper Fig. 11b).
CSX = AddDump(CSX, 'Hf', 'DumpType', 11, 'Frequency', [p.f_design]);
CSX = AddBox(CSX, 'Hf', 10, [-Lx, -Ly, -5*h], [Lx, Ly, 5*h]);

% NF2FF box kept in air, >=3 cells inside the MUR boundaries.
targets = [-45, 45, -45, 45, -15, 45];
[~, ix1] = min(abs(mesh.x - targets(1))); [~, ix2] = min(abs(mesh.x - targets(2)));
[~, iy1] = min(abs(mesh.y - targets(3))); [~, iy2] = min(abs(mesh.y - targets(4)));
[~, iz1] = min(abs(mesh.z - targets(5))); [~, iz2] = min(abs(mesh.z - targets(6)));
nfStart = [mesh.x(ix1), mesh.y(iy1), mesh.z(iz1)];
nfStop  = [mesh.x(ix2), mesh.y(iy2), mesh.z(iz2)];
[CSX, nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', nfStart, nfStop);

%% 8. Write, preview and run
CleanupSimPath(p.Sim_Path);
WriteOpenEMS(fullfile(p.Sim_Path, p.Sim_CSX), FDTD, CSX);
if p.show_geometry
    CSXGeomPlot(fullfile(p.Sim_Path, p.Sim_CSX));
end
if p.run_simulation
    RunOpenEMS(p.Sim_Path, p.Sim_CSX);
    % Multicore example (uncomment if your openEMS build supports it):
    % RunOpenEMS(p.Sim_Path, p.Sim_CSX, '--engine=multithreaded --numThreads=4');
else
    disp('run_simulation = false: XML written, solver skipped.');
    return;
end

%% 9. S-parameters, VSWR and feed impedance (1-4 GHz)
freq = linspace(1e9, 4e9, 601);
port = calcPort(port, p.Sim_Path, freq);
s11  = port.uf.ref ./ port.uf.inc;
S11dB = 20*log10(abs(s11));
VSWR  = (1 + abs(s11)) ./ (1 - abs(s11));
Zin   = port.uf.tot ./ port.if.tot;

[~, iRes] = min(abs(s11));      % magnitude-based resonance detection
f_res = freq(iRes);
disp(['Resonance: f = ' num2str(f_res/1e9, '%.4f') ' GHz, S11 = ' ...
      num2str(S11dB(iRes), '%.2f') ' dB, VSWR = ' ...
      num2str(VSWR(iRes), '%.3f')]);

figure('Name', 'Reflection coefficient S11');
plot(freq/1e6, S11dB, 'k-', 'LineWidth', 2);
grid on;
title('Reflection coefficient S_{11} (paper target: -25.81 dB @ 2.83 GHz)');
xlabel('frequency f / MHz');
ylabel('|S_{11}| / dB');

figure('Name', 'VSWR');
plot(freq/1e6, VSWR, 'b-', 'LineWidth', 2);
grid on;
title('VSWR at the 50 Ohm probe feed');
xlabel('frequency f / MHz');
ylabel('VSWR');

figure('Name', 'Feed point impedance');
plot(freq/1e6, real(Zin), 'k-', 'LineWidth', 2);
hold on;
plot(freq/1e6, imag(Zin), 'r--', 'LineWidth', 2);
grid on;
title('Feed point impedance');
xlabel('frequency f / MHz');
ylabel('impedance Z_{in} / Ohm');
legend('real', 'imag');

% Smith chart drawn manually (no toolbox / version dependency).
figure('Name', 'Smith chart');
hold on; axis equal; grid on;
th = linspace(0, 2*pi, 361);
plot(cos(th), sin(th), 'k-', 'LineWidth', 1.5);
for r = [0.2, 0.5, 1, 2, 5]
    plot(r/(r+1) + cos(th)/(r+1), sin(th)/(r+1), 'k:');
end
for x = [0.5, 1, 2, -0.5, -1, -2]
    plot(1 + cos(th)/x, 1/x + sin(th)/x, 'k:');
end
plot(real(s11), imag(s11), 'r-', 'LineWidth', 1.5);
plot(real(s11(iRes)), imag(s11(iRes)), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
title('Input reflection coefficient (o = resonance)');
xlabel('real'); ylabel('imag');
drawnow;

%% 10. Far field at resonance: directivity, gain, efficiency
disp('Calculating far field at phi = [0 90] deg ...');
nf2ff = CalcNF2FF(nf2ff, p.Sim_Path, f_res, ...
    [-180:2:180]*pi/180, [0, 90]*pi/180);

Pacc = port.P_inc(:) .* (1 - abs(s11(:)).^2);   % accepted power
eta_rad = nf2ff.Prad / Pacc(iRes);              % radiation efficiency
gain_lin = eta_rad * nf2ff.Dmax;                % G = eta_rad * D
disp(['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp(['directivity: Dmax = ' num2str(10*log10(nf2ff.Dmax), '%.2f') ...
      ' dBi (paper: 6.6 dB)']);
disp(['gain: G = ' num2str(10*log10(gain_lin), '%.2f') ...
      ' dB (paper: 7.3 dB)']);
disp(['efficiency: eta = ' num2str(100*eta_rad, '%.1f') ...
      ' % (paper: 96 %)']);

figure('Name', 'Far-field cuts (polar)');
polarFF(nf2ff, 'xaxis', 'theta', 'param', [1, 2], 'normalize', 1);

figure('Name', 'Far-field cuts (dB)');
plotFFdB(nf2ff, 'xaxis', 'theta', 'param', [1, 2]);
drawnow;

disp('Calculating 3D pattern (ParaView: tmp_CRM_Star_Patch/3D_Pattern.vtk) ...');
thetaRange = 0:2:180;
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, p.Sim_Path, f_res, thetaRange*pi/180, ...
    phiRange*pi/180, 'Verbose', 1, 'Outfile', '3D_Pattern.h5');
figure('Name', '3D far-field pattern');
plotFF3D(nf2ff, 'logscale', -20);

E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK(fullfile(p.Sim_Path, '3D_Pattern.vtk'), E_far_normalized, ...
    thetaRange, phiRange, 'scale', 1e-3);

save(fullfile(p.Sim_Path, 'crm_star_results.mat'), ...
    'freq', 's11', 'S11dB', 'VSWR', 'Zin', 'f_res', ...
    'eta_rad', 'gain_lin');
disp('Results saved to tmp_CRM_Star_Patch/crm_star_results.mat');
