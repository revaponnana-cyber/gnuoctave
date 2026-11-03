%% Simple Patch Antenna Tutorial
%% -----------------------------
%% This tutorial simulates a microstrip patch antenna on a low-loss dielectric
%% substrate, covering the full workflow from geometry setup to S-parameter
%% extraction and far-field pattern analysis. It introduces the key openEMS
%% building blocks — CSXCAD geometry, lumped port excitation, automatic meshing,
%% and the near-field-to-far-field transformation — in a single self-contained
%% example.
%
% Tested with
%  - Octave 11.3
%  - openEMS v0.37
%
% (C) 2010-2026 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% Setup the Simulation
%% --------------------
%% Define the physical dimensions and material properties of the patch,
%% substrate, and ground plane. The patch ``width`` controls the resonant
%% frequency — a narrower patch resonates higher — while the substrate
%% permittivity ``epsR`` and thickness together set the effective wavelength
%% inside the dielectric. The ``SimBox`` must be large enough that the
%% absorbing boundaries sit well away from the antenna near-field region.
physical_constants;
unit = 1e-3; % all length in mm

% patch width in x-direction
patch.width  = 44; % resonant length
% patch length in y-direction
patch.length = 55;

%substrate setup
% FR4 substrate parameters
substrate.epsR = 4.4;
substrate.tanD = 0.025;
substrate.width = 50;
substrate.length = 50;
substrate.thickness = 1.6;
substrate.cells = 5;

% Dielectric loss
substrate.kappa = 2*pi*2.82e9*EPS0*substrate.epsR*substrate.tanD;

%setup feeding
% Verified clear of all slots (nearest slot edge >1.1 mm away) after the
% star rotation. Tune along x in [-8,-4] for the deepest S11 null.
feed.pos = -6; %feeding position in x-direction
feed.R = 50;     %feed resistance

% size of the simulation box
SimBox = [200 200 150];

%% Setup FDTD Parameter & Excitation Function
%% ------------------------------------------
%% Configure the broadband Gaussian pulse and the first-order Mur absorbing
%% boundary conditions. The centre frequency ``f0`` and 20 dB bandwidth ``fc``
%% should bracket the expected patch resonance so that the full S-parameter
%% response is captured in a single simulation run. Mur boundaries are chosen
%% here because the antenna sits well inside the box; switch to PML if you
%% need stronger absorption for a more compact simulation domain.
f0 = 2.5e9; % center frequency
fc = 1.5e9; % 20 dB corner frequency
% NrTS is only a cap: the run stops early once the energy has decayed by
% EndCriteria. 1e-4 (-40 dB) is the standard S11-grade threshold.
FDTD = InitFDTD('NrTS', 200000, 'EndCriteria', 1e-4);
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup CSXCAD Geometry & Mesh
%% ----------------------------
%% Build the complete antenna structure — patch metal, dielectric substrate,
%% ground plane, lumped feed port, and the NF2FF integration surface — inside
%% a single CSXCAD container. The mesh is generated with ``DetectEdges`` and
%% ``SmoothMesh``: patch edges are resolved at lambda/50 to capture the patch
%% current distribution accurately, while the surrounding free-space region
%% uses the coarser lambda/20 spacing to keep the total cell count manageable.
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/3 SimBox(3)*2/3];


% Create Circular Patch
CSX = AddMetal(CSX,'patch');

patch.radius = 14.5;

% Create Air material for star-shaped slots
CSX = AddMaterial(CSX, 'Air');
CSX = SetMaterialProperty(CSX, 'Air', 'Epsilon', 1, 'Mue', 1);

% Create 100 unique points around the circle
theta = (0:99)*2*pi/100;

points = [patch.radius*cos(theta);
          patch.radius*sin(theta)];

CSX = AddPolygon(CSX,'patch',10,'z',substrate.thickness,points);
% Create first six-pointed star outline slot
% Regular hexagram (inner = outer/sqrt(3), tips pointing up) matching the
% paper shape. The 30-deg rotation also moves every slot off the feed point
% at (-6,0), which used to sit exactly on a star tip vertex.
star.outerRadius = 8.0;
star.innerRadius = 8.0/sqrt(3);
star.slotWidth = 0.5;
star.rot = pi/6;

% Generate the 12 alternating outer and inner star points
starPoints = zeros(2,12);

for k = 1:12
    angle = (k-1)*pi/6 + star.rot;

    if mod(k,2) == 1
        r = star.outerRadius;
    else
        r = star.innerRadius;
    end

    starPoints(:,k) = [r*cos(angle); r*sin(angle)];
end

% Cut a thin air slot along each star edge
for k = 1:12
    next = mod(k,12) + 1;

    p1 = starPoints(:,k);
    p2 = starPoints(:,next);

    edge = p2 - p1;
    edgeLength = norm(edge);
    normal = [-edge(2); edge(1)] / edgeLength;

    offset = normal * star.slotWidth/2;

    slotPoints = [p1+offset, p2+offset, ...
                  p2-offset, p1-offset];

    CSX = AddPolygon(CSX,'Air',20,'z', ...
                     substrate.thickness,slotPoints);
end
% Create second, smaller six-pointed star outline slot
star2.outerRadius = 6.0;
star2.innerRadius = 6.0/sqrt(3);
star2.slotWidth = 0.5;
star2.rot = pi/6;

star2Points = zeros(2,12);

for k = 1:12
    angle = (k-1)*pi/6 + star2.rot;

    if mod(k,2) == 1
        r = star2.outerRadius;
    else
        r = star2.innerRadius;
    end

    star2Points(:,k) = [r*cos(angle); r*sin(angle)];
end

for k = 1:12
    next = mod(k,12) + 1;

    p1 = star2Points(:,k);
    p2 = star2Points(:,next);

    edge = p2 - p1;
    edgeLength = norm(edge);
    normal = [-edge(2); edge(1)] / edgeLength;

    offset = normal * star2.slotWidth/2;

    slotPoints = [p1+offset, p2+offset, ...
                  p2-offset, p1-offset];

    CSX = AddPolygon(CSX,'Air',20,'z', ...
                     substrate.thickness,slotPoints);
end

% Create third, smallest six-pointed star outline slot
star3.outerRadius = 4.0;
star3.innerRadius = 4.0/sqrt(3);
star3.slotWidth = 0.5;
star3.rot = pi/6;

star3Points = zeros(2,12);

for k = 1:12
    angle = (k-1)*pi/6 + star3.rot;

    if mod(k,2) == 1
        r = star3.outerRadius;
    else
        r = star3.innerRadius;
    end

    star3Points(:,k) = [r*cos(angle); r*sin(angle)];
end

for k = 1:12
    next = mod(k,12) + 1;

    p1 = star3Points(:,k);
    p2 = star3Points(:,next);

    edge = p2 - p1;
    edgeLength = norm(edge);
    normal = [-edge(2); edge(1)] / edgeLength;

    offset = normal * star3.slotWidth/2;

    slotPoints = [p1+offset, p2+offset, ...
                  p2-offset, p1-offset];

    CSX = AddPolygon(CSX,'Air',20,'z', ...
                     substrate.thickness,slotPoints);
end
% Create thin circular outline around the three stars
ring.radius = 10.5;
ring.slotWidth = 0.5;
ring.segments = 120;

for k = 1:ring.segments
    angle1 = (k-1)*2*pi/ring.segments;
    angle2 = k*2*pi/ring.segments;
        % Leave small gaps at the top and bottom of the ring
    angleMid = (angle1 + angle2)/2;
    gapHalfAngle = (1.0/2)/ring.radius;

    topDifference = abs(atan2( ...
        sin(angleMid - pi/2), cos(angleMid - pi/2)));

    bottomDifference = abs(atan2( ...
        sin(angleMid - 3*pi/2), cos(angleMid - 3*pi/2)));

    if topDifference < gapHalfAngle || ...
       bottomDifference < gapHalfAngle
        continue;
    end

    outerRadius = ring.radius + ring.slotWidth/2;
    innerRadius = ring.radius - ring.slotWidth/2;

    ringPoints = [
        outerRadius*cos(angle1), outerRadius*cos(angle2), ...
        innerRadius*cos(angle2), innerRadius*cos(angle1);
        outerRadius*sin(angle1), outerRadius*sin(angle2), ...
        innerRadius*sin(angle2), innerRadius*sin(angle1)
    ];

    CSX = AddPolygon(CSX,'Air',20,'z', ...
                     substrate.thickness,ringPoints);
end

% ==========================================
% CRM DESIGN - FIRST UNIT
% ==========================================

% CRM centre position
% Four CRM corner positions
crmPositions = [
    -18,  18;   % Top-left
     18,  18;   % Top-right
    -18, -18;   % Bottom-left
     18, -18    % Bottom-right
];

for crmIndex = 1:4
    crm.centerX = crmPositions(crmIndex,1);
    crm.centerY = crmPositions(crmIndex,2);

% Dimensions from research paper Table 2
crm.M = 0.2;
crm.T = 5.0;
crm.U = 5.95;

% Ring drawing resolution
crm.segments = 120;

% Draw the outer circular ring
for k = 1:crm.segments

    angle1 = (k-1)*2*pi/crm.segments;
    angle2 = k*2*pi/crm.segments;

    rOuter = crm.U;
    rInner = crm.U - crm.M;

    ringPoints = [
        rOuter*cos(angle1), rOuter*cos(angle2), ...
        rInner*cos(angle2), rInner*cos(angle1);
        rOuter*sin(angle1), rOuter*sin(angle2), ...
        rInner*sin(angle2), rInner*sin(angle1)
    ];

    ringPoints(1,:) = ringPoints(1,:) + crm.centerX;
    ringPoints(2,:) = ringPoints(2,:) + crm.centerY;

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,ringPoints);
end

% Draw the inner circular ring
% Draw the inner circular ring with four openings

for k = 1:crm.segments

    angle1 = (k-1)*2*pi/crm.segments;
    angle2 = k*2*pi/crm.segments;
    angleMid = (angle1 + angle2)/2;

    % Leave openings at 12, 3, 6 and 9 o'clock
    gapHalfAngle = 0.25/crm.T;

    skipSegment = false;

    for j = 0:3
        gapAngle = j*pi/2;

        angleDifference = abs(atan2( ...
            sin(angleMid-gapAngle), ...
            cos(angleMid-gapAngle)));

        if angleDifference < gapHalfAngle
            skipSegment = true;
        end
    end

    if skipSegment
        continue;
    end

    rOuter = crm.T;
    rInner = crm.T - crm.M;

    ringPoints = [
        rOuter*cos(angle1), rOuter*cos(angle2), ...
        rInner*cos(angle2), rInner*cos(angle1);
        rOuter*sin(angle1), rOuter*sin(angle2), ...
        rInner*sin(angle2), rInner*sin(angle1)
    ];

    ringPoints(1,:) = ringPoints(1,:) + crm.centerX;
    ringPoints(2,:) = ringPoints(2,:) + crm.centerY;

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,ringPoints);
end
% Four large circular pads at cardinal positions

crm.dumbbellRadius = 1.5;  % R from paper
crm.dumbbellDistance = 2.2; % Estimated centre position

for k = 1:4

    angle = k*pi/2;

    padX = crm.centerX + crm.dumbbellDistance*cos(angle);
    padY = crm.centerY + crm.dumbbellDistance*sin(angle);

    padPoints = zeros(2,crm.segments);

    for n = 1:crm.segments

        a = (n-1)*2*pi/crm.segments;

        padPoints(:,n) = [
            padX + crm.dumbbellRadius*cos(a);
            padY + crm.dumbbellRadius*sin(a)
        ];

    end

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,padPoints);

end
% Two parallel connecting strips for each large circular pad

crm.dumbbellStripWidth = 0.3;  % Estimated strip width
crm.dumbbellGap = 0.3;         % Estimated gap between strips

for k = 0:3

    angle = k*pi/2;

    % Start at the outside edge of the large circular pad
   rStart = crm.dumbbellDistance + crm.dumbbellRadius - 0.4;

    % End at the inner-ring opening
    rEnd = crm.T - crm.M;

    % Draw two strips, one on each side of the centre line
    for side = [-1 1]

        offset = side * ...
            (crm.dumbbellGap/2 + crm.dumbbellStripWidth/2);

        halfWidth = crm.dumbbellStripWidth/2;

        localPoints = [
            offset-halfWidth, offset+halfWidth, ...
            offset+halfWidth, offset-halfWidth;
            rStart, rStart, rEnd, rEnd
        ];

        rotation = [
            cos(angle), -sin(angle);
            sin(angle),  cos(angle)
        ];

        stripPoints = rotation*localPoints;

        stripPoints(1,:) = stripPoints(1,:) + crm.centerX;
        stripPoints(2,:) = stripPoints(2,:) + crm.centerY;

        CSX = AddPolygon(CSX,'patch',20,'z', ...
            substrate.thickness,stripPoints);

    end
end
% Four diagonal lollipop elements

crm.lollipopRadius = 0.8;  % N: circular head radius
crm.lollipopDistance = 7.0;
crm.lollipopTailWidth = 0.3;
crm.lollipopTailEnd = 5.3;

for k = 1:4

    angle = pi/4 + (k-1)*pi/2;

    % Centre of circular head
    headX = crm.centerX + crm.lollipopDistance*cos(angle);
    headY = crm.centerY + crm.lollipopDistance*sin(angle);

    % Draw circular head
    headPoints = zeros(2,crm.segments);

    for n = 1:crm.segments
        a = (n-1)*2*pi/crm.segments;

        headPoints(:,n) = [
            headX + crm.lollipopRadius*cos(a);
            headY + crm.lollipopRadius*sin(a)
        ];
    end

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,headPoints);

    % Draw narrow tail extending inward
    rStart = crm.lollipopDistance;
    rEnd = crm.lollipopTailEnd;
    halfWidth = crm.lollipopTailWidth/2;

    localPoints = [
        -halfWidth, halfWidth, halfWidth, -halfWidth;
        rEnd, rEnd, rStart, rStart
    ];

    rotation = [
        cos(angle), -sin(angle);
        sin(angle),  cos(angle)
    ];

    tailPoints = rotation*localPoints;

    tailPoints(1,:) = tailPoints(1,:) + crm.centerX;
    tailPoints(2,:) = tailPoints(2,:) + crm.centerY;

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,tailPoints);
end
% Four inward-facing metallic strips

crm.stubWidth = 0.3;
crm.stubOuterRadius = crm.U;
crm.stubInnerRadius = crm.T + 0.3;

for k = 1:4

    angle = (k-1)*pi/2;

    r1 = crm.stubOuterRadius;
    r2 = crm.stubInnerRadius;
    halfWidth = crm.stubWidth/2;

    localPoints = [
        -halfWidth, halfWidth, halfWidth, -halfWidth;
        r2, r2, r1, r1
    ];

    rotation = [
        cos(angle), -sin(angle);
        sin(angle),  cos(angle)
    ];

    stubPoints = rotation*localPoints;

    stubPoints(1,:) = stubPoints(1,:) + crm.centerX;
    stubPoints(2,:) = stubPoints(2,:) + crm.centerY;

    CSX = AddPolygon(CSX,'patch',20,'z', ...
        substrate.thickness,stubPoints);
end
end % End of four-CRM placement loop

% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );
start = [-substrate.width/2 -substrate.length/2 0];
stop  = [ substrate.width/2  substrate.length/2 substrate.thickness];
CSX = AddBox( CSX, 'substrate', 0, start, stop );

% add extra cells to discretize the substrate thickness
% Refine mesh through substrate and add air-region mesh above patch
mesh.z = [mesh.z, linspace(0, substrate.thickness, 6), ...
          2, 5, 10, 15, 20, 30, 40, 60, 80];

% air-region lines: below the ground (for the NF2FF bottom face) and around
% the antenna in x/y (coarse air cells plus faces for the NF2FF box)
mesh.z = [mesh.z, -40, -30, -20, -12, -6, -2];
mesh.x = [mesh.x, -85, -70, -55, -40, -30, 30, 40, 55, 70, 85];
mesh.y = [mesh.y, -85, -70, -55, -40, -30, 30, 40, 55, 70, 85];

mesh.z = uniquetol(mesh.z, 1e-10);

% Create Ground same size as substrate
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start(3)=0;
stop(3) =0;
CSX = AddBox(CSX,'gnd',10,start,stop);

% Cross-shaped defected ground structure (paper: J=20 x K=1.2 arms).
% Improves the impedance match and radiation efficiency; set DGS_enabled
% to false for a solid ground plane (ablation case).
DGS_enabled = true;
if DGS_enabled
    dgs.arm = 20; dgs.w = 1.2; dgs.cx = 0; dgs.cy = 14;
    CSX = AddBox(CSX,'Air',20, ...
        [dgs.cx-dgs.arm/2 dgs.cy-dgs.w/2 -0.5], ...
        [dgs.cx+dgs.arm/2 dgs.cy+dgs.w/2 0.2]);
    CSX = AddBox(CSX,'Air',20, ...
        [dgs.cx-dgs.w/2 dgs.cy-dgs.arm/2 -0.5], ...
        [dgs.cx+dgs.w/2 dgs.cy+dgs.arm/2 0.2]);
end

% Apply the Excitation & Resist as a Current Source
start = [feed.pos 0 0];
stop  = [feed.pos 0 substrate.thickness];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

% Finalize the Mesh
% detect all edges except of the patch
mesh = DetectEdges(CSX, mesh,'ExcludeProperty','patch');
% detect and set a special 2D metal edge mesh for the patch
mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', 1);
% Remove duplicate mesh coordinates
% Remove duplicate and nearly duplicate mesh coordinates
mesh.x = uniquetol(mesh.x, 1e-10);
mesh.y = uniquetol(mesh.y, 1e-10);
mesh.z = uniquetol(mesh.z, 1e-10);
% generate a smooth mesh with max. cell size: lambda_min / 20
% (fills the large air gaps with graded cells so waves are resolved and the
% NF2FF box sits several cells inside the MUR boundaries)
mesh = SmoothMesh(mesh, c0/(f0+fc)/unit/20);

disp('Mesh diagnostics:');
disp(['Number of X mesh lines: ', num2str(length(mesh.x))]);
disp(['Number of Y mesh lines: ', num2str(length(mesh.y))]);
disp(['Number of Z mesh lines: ', num2str(length(mesh.z))]);

disp(['Smallest X spacing: ', num2str(min(diff(sort(mesh.x))))]);
disp(['Smallest Y spacing: ', num2str(min(diff(sort(mesh.y))))]);
disp(['Smallest Z spacing: ', num2str(min(diff(sort(mesh.z))))]);

disp('Duplicate or nearly duplicate X mesh coordinates:');
dx = diff(sort(mesh.x));
disp(find(dx < 1e-10));

disp('Duplicate or nearly duplicate Y mesh coordinates:');
dy = diff(sort(mesh.y));
disp(find(dy < 1e-10));

CSX = DefineRectGrid(CSX, unit, mesh);

CSX = AddDump(CSX,'Hf', 'DumpType', 11, 'Frequency',[2.82e9]); % surface-current view at resonance
CSX = AddBox(CSX,'Hf',10,[-substrate.width -substrate.length -10*substrate.thickness],[substrate.width +substrate.length 10*substrate.thickness]); %assign box

% add a nf2ff calc box; size is 3 cells away from MUR boundary condition
% Create NF2FF box around the antenna
% Create NF2FF box around the antenna with air clearance
% NF2FF faces snapped to the nearest mesh lines: the box encloses the whole
% antenna (substrate/CRM extend to ~25 mm) and stays in air, well clear of
% the MUR boundaries at x/y=+-100, z=-50/+100.
[~, ix1] = min(abs(mesh.x - (-40))); [~, ix2] = min(abs(mesh.x - 40));
[~, iy1] = min(abs(mesh.y - (-40))); [~, iy2] = min(abs(mesh.y - 40));
[~, iz1] = min(abs(mesh.z - (-12))); [~, iz2] = min(abs(mesh.z - 40));
start = [mesh.x(ix1) mesh.y(iy1) mesh.z(iz1)];
stop  = [mesh.x(ix2) mesh.y(iy2) mesh.z(iz2)];

[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);
disp('NF2FF start coordinates:');
disp(start);

disp('NF2FF stop coordinates:');
disp(stop);

disp('Z mesh coordinates:');
disp(mesh.z);

%% Prepare and Run Simulation
%% --------------------------
%% Write the geometry and solver settings to an XML file, optionally preview
%% the structure in AppCSXCAD, then launch the openEMS solver. Inspecting the
%% geometry plot before running is strongly recommended: it lets you verify
%% patch, substrate, and port placement immediately, without waiting for the
%% solver to finish.
Sim_Path = 'tmp_Patch_Ant';
Sim_CSX = 'patch_ant.xml';

% create an empty working directory
CleanupSimPath(Sim_Path);

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% set RUN_OPENEMS=false to only write/check geometry without simulating
RUN_OPENEMS = true;
if ~RUN_OPENEMS
    disp('RUN_OPENEMS=false: geometry written, solver skipped.');
    return;
end


% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% Postprocessing & Plots
%% ----------------------
%% Read back the time-domain port voltages and currents recorded during the
%% simulation and transform them to the frequency domain. ``calcPort`` computes
%% the incident and reflected wave amplitudes at each frequency point; these
%% are the basis for the impedance, S11, and radiation-efficiency calculations
%% that follow.
freq = linspace(1e9, 4e9, 601);
port = calcPort(port, Sim_Path, freq);

%% Smith Chart Port Reflection
%% ---------------------------
%% Plot the reflection coefficient on a Smith chart and as an S11 magnitude
%% curve to characterise the antenna input match. A deep null in ``|S11|`` at the
%% design frequency confirms resonance; the -10 dB bandwidth indicates the
%% usable frequency range over which the antenna is well matched to the
%% 50 Ohm feed.
try
    plotRefl(port, 'threshold', -10);
catch
    try
        plotRefl(port, freq, 'threshold', -10);
    catch
        disp('plotRefl not available in this openEMS version; skipping Smith chart.');
    end
end
title( 'reflection coefficient' );

% plot feed point impedance
Zin = port.uf.tot ./ port.if.tot;
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
s11 = port.uf.ref ./ port.uf.inc;
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

% plot VSWR (paper Fig. 7)
VSWR = (1+abs(s11))./(1-abs(s11));
figure
plot( freq/1e6, VSWR, 'b-', 'Linewidth', 2 );
grid on
title( 'VSWR at the 50 Ohm feed' );
xlabel( 'frequency f / MHz' );
ylabel( 'VSWR' );

drawnow

%% NFFF Plots
%% ----------
%% Invoke the near-field-to-far-field transformation to obtain the radiation
%% pattern at the resonant frequency identified from the S11 minimum. Sampling
%% at phi = 0 and 90 degrees reveals the two principal-plane cuts; the
%% subsequent full-sphere calculation yields maximum directivity and radiation
%% efficiency, showing how effectively the antenna converts accepted power into
%% useful far-field radiation.
%find resonance frequency from s11 (deepest |S11| null)
[~, f_res_ind] = min(abs(s11));
f_res = freq(f_res_ind);

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, [-180:2:180]*pi/180, [0 90]*pi/180);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
P_acc = port.P_acc(f_res_ind); % accepted power from calcPort
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./P_acc) ' %']);
gain_lin = (nf2ff.Prad./P_acc) * nf2ff.Dmax; % G = eta_rad * D
disp( ['gain: G = ' num2str(10*log10(gain_lin)) ' dB (paper: 7.3 dB)'] );

% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',[1 2],'normalize',1)

% log-scale directivity plot
figure
plotFFdB(nf2ff,'xaxis','theta','param',[1 2])
% conventional plot approach
% plot( nf2ff.theta*180/pi, 20*log10(nf2ff.E_norm{1}/max(nf2ff.E_norm{1}(:)))+10*log10(nf2ff.Dmax));

drawnow

% Show 3D pattern
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);


E_far_normalized = nf2ff.E_norm{1} / max(nf2ff.E_norm{1}(:)) * nf2ff.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);

% save port/nf2ff structs for PlotStarSlotResults.m (standalone viewer)
save([Sim_Path '/results_workspace.mat'], 'port', 'nf2ff', 'freq', 's11', 'f_res', 'f_res_ind');
disp('Saved results_workspace.mat for the standalone results viewer.');
