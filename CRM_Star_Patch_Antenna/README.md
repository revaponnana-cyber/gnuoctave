# CRM-inspired star-shaped microstrip patch antenna (MATLAB / openEMS)

Full-wave model of the paper design: a circular FR-4 patch (50 x 50 x 1.6 mm³,
εr = 4.4, tanδ = 0.025) with triple star slots and an outer ring slot, four
corner circular-ring-metamaterial (CRM) unit cells, a cross-shaped defected
ground structure (Cr-DGS), and a 50 Ω coaxial-probe feed. Design frequency:
**2.82 GHz**, analysed over **1-4 GHz**.

Paper targets (simulated): S11 = **-25.81 dB @ 2.83 GHz**, gain **7.3 dB**,
directivity **6.6 dB**, radiation efficiency **96 %**.

> Note: the code below is an independent implementation built from the paper's
> tables and figures; the paper itself (© Springer) is not redistributed here.

## Files

| File | Needs openEMS? | Purpose |
|---|---|---|
| `CRM_Star_Patch_Antenna.m` | yes | Main script: geometry, mesh, FDTD run, S11/VSWR/Zin/Smith, NF2FF gain/directivity/efficiency, 3D pattern + VTK |
| `antenna_params.m` | no | All dimensions (paper Tables 1-2) + FDTD/run control in one struct |
| `hexagram_outline.m`, `arc_band.m`, `circle_points.m` | no | 2D outline generators (star centre-lines, ring sectors, disks) |
| `staircase_boxes.m`, `add_crm_cell.m` | yes | CSXCAD builders: slot staircases, one CRM unit cell |
| `snap_coords.m`, `snap_polygon.m` | no | Mesh-grid quantization helpers |
| `preview_layout.m` | **no** | 2D top/bottom layout plot - run first, anywhere |
| `circular_patch_theory.m` | **no** | Paper Eqs. (1)-(2): radius design curve + radius-sweep table |

## Quick start

```matlab
cd CRM_Star_Patch_Antenna
preview_layout            % 1. check the geometry (no openEMS needed)
circular_patch_theory     % 2. analytical radius check (prints ~14.4 mm)
CRM_Star_Patch_Antenna    % 3. full-wave simulation (needs openEMS)
```

Requirements for step 3: MATLAB R2016b+ (or GNU Octave) with the openEMS +
CSXCAD MATLAB/Octave interface installed and on the path, plus the `openEMS`
solver binary. Set `p.run_simulation = false` in `antenna_params.m` to only
write `tmp_CRM_Star_Patch/crm_star_patch.xml` and preview the geometry.

## Ablation studies (paper Tables 3-4)

```matlab
p.crm.enabled = false;   % "without CRM" case  -> gain ~2.8 dB, eta ~85 %
p.dgs.enabled = false;   % solid ground plane (no cross DGS)
```

## Dimension mapping and assumptions

Table 1 (mm): A=14.5, B=4.6, C=0.75, D=4.04, E=1.84, F=1.24, G=2.43, H=6,
I=6.66, J=20, K=1.2, L=1. Table 2 (mm): M=0.2, N=0.8, O=0.95, P=0.44, Q=0.3,
R=1.5, S=0.95, T=5, U=5.95.

| Param | Used as | Confidence |
|---|---|---|
| A | patch radius 14.5 | confirmed (paper text, §4.1.2) |
| B, D, E | star tip-to-notch x2 → tip radii [9.2, 8.08, 3.68] | figure-derived |
| C | ring-slot width 0.75; also outer/mid star width (keeps nested slots separate) | figure-derived |
| F | inner-star slot width 1.24 | figure-derived |
| G | ring-slot opening width 2.43 at bottom | assumed |
| H, I | cross-check only (predicted ≈ 5.9/6.6 vs 6/6.66) | - |
| J / K | cross arms 20 x 1.2, centred (0, 14) | figure scale |
| L | pin radius 1 (2 mm dot) | figure scale |
| T, U / M / Q | ring radii 5/5.95, trace 0.2→0.25, splits 0.3→0.25 (mesh-grid quantized, ≤0.13 mm shift) | table + mesh |
| R / N | inner/outer disk **diameters** 1.5/0.8 | figure scale |
| S / P | stub length 0.95 / width 0.44→0.5 | figure-derived |
| O | ring pitch U−T = 0.95 (descriptive) | exact |
| CRM centres | (±16.8, ±16.8) | assumed (Fig. 2a) |
| feed | (0, −5.5): Fig. 4b shows ≈ −3.5, shifted to clear the inner-star slot | assumed + tuned |

Modelling choices: slots are cut by higher-priority air boxes over the patch /
ground sheets; the coax probe is a 50 Ω vertical lumped port (standard openEMS
practice); diagonal edges use grid-snapped staircase boxes (0.25 mm) to bound
the mesh to roughly 150 x 150 x 25 cells. The "two semi-circular side cuts"
mentioned once in the text are not visible in any figure and are not modelled.

## Expected mesh / runtime

~0.4-1.0 M cells, dt ≈ 0.48 ps, ≤ 50 000 steps with early stop at 1e-4 energy
decay. Typically tens of minutes with a multithreaded openEMS build (see the
`RunOpenEMS` example in the main script). Reduce `p.SimBox` or raise
`p.EndCriteria` for faster draft runs.

## Paper quirks to be aware of

- Reported gain (7.3 dB) exceeds reported directivity (6.6 dB), which violates
  G = η·D ≤ D; likely a reporting/labelling slip. This script derives both
  from one NF2FF solution, so gain ≤ directivity always holds.
- "VSWR of nearest 2 dB" (§4.1): VSWR is unitless (and ≈ 2 would contradict
  S11 = −25.8 dB → VSWR ≈ 1.1). The script plots true VSWR.
- Eq. (2) uses a +0.5 exponent where the classic microstrip form has −0.5;
  `circular_patch_theory.m` implements the paper literally and still lands at
  ≈ 14.4 mm ≈ A, confirming the patch size.
