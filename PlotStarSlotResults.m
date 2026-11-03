%% PlotStarSlotResults -- standalone results viewer (no FDTD re-run)
%% ------------------------------------------------------------------
%% Replots every result (S11, VSWR, impedance, Smith chart, gain /
%% directivity / efficiency vs frequency, radiation patterns) from an
%% existing run in tmp_Patch_Ant, and prints a summary table with the
%% paper targets for comparison. Figures are also saved as PNGs plus a
%% results_summary.txt, so you can inspect the current results anytime.
%%
%% Workflow in GNU Octave (with openEMS on the path):
%%   1. Star_Patch_Antenna_StarSlots   % full run once (slow; saves data)
%%   2. PlotStarSlotResults            % re-view results anytime (fast)
%%
%% Paper targets (simulated): S11 = -25.81 dB @ 2.83 GHz, gain 7.3 dB,
%% directivity 6.6 dB, efficiency 96 %.

close all
clear
clc

Sim_Path  = 'tmp_Patch_Ant';
WS_FILE   = [Sim_Path '/results_workspace.mat'];
SAVE_FIGS = true;   % also write PNGs + summary txt into results_figs/

if ~exist(WS_FILE, 'file')
    error(['No ' WS_FILE ' found. Run Star_Patch_Antenna_StarSlots once first.']);
end
load(WS_FILE);      % port + nf2ff structs saved by the finished run

%% 1. Port results over 1-4 GHz (re-read from recorded time data)
freq = linspace(1e9, 4e9, 601);
port = calcPort(port, Sim_Path, freq);
s11   = port.uf.ref ./ port.uf.inc;
S11dB = 20*log10(abs(s11));
VSWR  = (1 + abs(s11)) ./ (1 - abs(s11));
Zin   = port.uf.tot ./ port.if.tot;

[~, iRes] = min(abs(s11));
f_res = freq(iRes);

% -10 dB bandwidth around the resonance
below10 = S11dB < -10;
i1 = iRes; while i1 > 1 && below10(i1-1), i1 = i1 - 1; end
i2 = iRes; while i2 < length(freq) && below10(i2+1), i2 = i2 + 1; end
if below10(iRes)
    BW = freq(i2) - freq(i1);
else
    BW = 0;
end

figdir = [Sim_Path '/results_figs'];
if SAVE_FIGS && ~exist(figdir, 'dir')
    mkdir(figdir);
end

%% 2. Plots: S11, VSWR, impedance, Smith chart
figure('Name', 'S11');
plot(freq/1e6, S11dB, 'k-', 'LineWidth', 2);
grid on; hold on;
plot([freq(1) freq(end)]/1e6, [-10 -10], 'r--');
plot(f_res/1e6, S11dB(iRes), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
title('Reflection coefficient S_{11} (paper: -25.81 dB @ 2.83 GHz)');
xlabel('frequency f / MHz'); ylabel('|S_{11}| / dB');
if SAVE_FIGS, save_png(figdir, 'S11.png'); end

figure('Name', 'VSWR');
plot(freq/1e6, VSWR, 'b-', 'LineWidth', 2);
grid on;
title('VSWR at the 50 Ohm feed');
xlabel('frequency f / MHz'); ylabel('VSWR');
if SAVE_FIGS, save_png(figdir, 'VSWR.png'); end

figure('Name', 'Impedance');
plot(freq/1e6, real(Zin), 'k-', 'LineWidth', 2);
hold on; grid on;
plot(freq/1e6, imag(Zin), 'r--', 'LineWidth', 2);
title('Feed point impedance');
xlabel('frequency f / MHz'); ylabel('impedance Z_{in} / Ohm');
legend('real', 'imag');
if SAVE_FIGS, save_png(figdir, 'impedance.png'); end

figure('Name', 'Smith');
hold on; axis equal; grid on;
thc = linspace(0, 2*pi, 361);
plot(cos(thc), sin(thc), 'k-', 'LineWidth', 1.5);
for r = [0.2 0.5 1 2 5]
    plot(r/(r+1) + cos(thc)/(r+1), sin(thc)/(r+1), 'k:');
end
for x = [0.5 1 2 -0.5 -1 -2]
    plot(1 + cos(thc)/x, 1/x + sin(thc)/x, 'k:');
end
plot(real(s11), imag(s11), 'r-', 'LineWidth', 1.5);
plot(real(s11(iRes)), imag(s11(iRes)), 'bo', ...
    'MarkerFaceColor', 'b', 'MarkerSize', 8);
title('Input reflection coefficient (o = resonance)');
xlabel('real'); ylabel('imag');
if SAVE_FIGS, save_png(figdir, 'smith.png'); end
drawnow;

%% 3. Gain / directivity / efficiency vs frequency (NF2FF sweep, no FDTD)
disp('Sweeping NF2FF over frequency (this reads the saved fields)...');
fsw   = linspace(2.0e9, 3.6e9, 17);
th_sw = (0:5:180)*pi/180;
ph_sw = (0:15:345)*pi/180;
Gsw = zeros(size(fsw)); Dsw = Gsw; Esw = Gsw;
Pacc_fine = port.P_inc(:) .* (1 - abs(s11(:)).^2);
for k = 1:length(fsw)
    disp(['  NF2FF ' num2str(k) '/' num2str(length(fsw)) ' @ ' ...
          num2str(fsw(k)/1e9, '%.3f') ' GHz']);
    nk = CalcNF2FF(nf2ff, Sim_Path, fsw(k), th_sw, ph_sw);
    Pacc_k = interp1(freq, Pacc_fine, fsw(k));
    eta_k = nk.Prad / Pacc_k;
    Gsw(k) = eta_k * nk.Dmax;
    Dsw(k) = nk.Dmax;
    Esw(k) = eta_k;
end

figure('Name', 'GainVsFreq');
subplot(3,1,1);
plot(fsw/1e9, 10*log10(Gsw), 'k-o', 'LineWidth', 1.5);
grid on; title('Peak gain vs frequency (paper: 7.3 dB)');
ylabel('gain / dB');
subplot(3,1,2);
plot(fsw/1e9, 10*log10(Dsw), 'b-o', 'LineWidth', 1.5);
grid on; title('Peak directivity vs frequency (paper: 6.6 dB)');
ylabel('directivity / dBi');
subplot(3,1,3);
plot(fsw/1e9, 100*Esw, 'r-o', 'LineWidth', 1.5);
grid on; title('Radiation efficiency vs frequency (paper: 96 %)');
ylabel('efficiency / %'); xlabel('frequency / GHz');
if SAVE_FIGS, save_png(figdir, 'gain_vs_freq.png'); end
drawnow;

%% 4. Radiation patterns at resonance
disp('Far-field cuts at resonance...');
nf2ff_res = CalcNF2FF(nf2ff, Sim_Path, f_res, ...
    [-180:2:180]*pi/180, [0 90]*pi/180);
figure('Name', 'FFpolar');
polarFF(nf2ff_res, 'xaxis', 'theta', 'param', [1 2], 'normalize', 1);
if SAVE_FIGS, save_png(figdir, 'pattern_polar.png'); end
figure('Name', 'FFdB');
plotFFdB(nf2ff_res, 'xaxis', 'theta', 'param', [1 2]);
if SAVE_FIGS, save_png(figdir, 'pattern_dB.png'); end
drawnow;

%% 5. Summary table: current results vs paper
G_res = interp1(fsw, Gsw, f_res);
D_res = nf2ff_res.Dmax;
eta_res = interp1(fsw, Esw, f_res);
[~, iG] = max(Gsw);
lines = {
    '=== Star-slot patch: current results ==='
    ['f_res      : ' num2str(f_res/1e9, '%.4f') ' GHz, S11 = ' ...
     num2str(S11dB(iRes), '%.2f') ' dB   (paper: -25.81 dB @ 2.83 GHz)']
    ['-10 dB BW  : ' num2str(BW/1e6, '%.1f') ' MHz']
    ['VSWR @ res : ' num2str(VSWR(iRes), '%.3f')]
    ['Zin @ res  : ' num2str(real(Zin(iRes)), '%.1f') ' + j' ...
     num2str(imag(Zin(iRes)), '%.1f') ' Ohm']
    ['Gain       : ' num2str(10*log10(G_res), '%.2f') ' dB @ f_res; peak ' ...
     num2str(10*log10(Gsw(iG)), '%.2f') ' dB @ ' ...
     num2str(fsw(iG)/1e9, '%.3f') ' GHz   (paper: 7.3 dB)']
    ['Directivity: ' num2str(10*log10(D_res), '%.2f') ...
     ' dBi @ f_res   (paper: 6.6 dB)']
    ['Efficiency : ' num2str(100*eta_res, '%.1f') ...
     ' % @ f_res   (paper: 96 %)']
};
for k = 1:length(lines), disp(lines{k}); end

if SAVE_FIGS
    fid = fopen([figdir '/results_summary.txt'], 'w');
    for k = 1:length(lines), fprintf(fid, '%s\n', lines{k}); end
    fclose(fid);
    save([Sim_Path '/results_table.mat'], 'freq', 'S11dB', 'VSWR', ...
        'Zin', 'f_res', 'fsw', 'Gsw', 'Dsw', 'Esw');
    disp(['Figures + summary saved in ' figdir]);
end

function save_png(figdir, name)
% save current figure as PNG; never abort the script if it fails
try
    print(gcf, [figdir '/' name], '-dpng');
catch
    disp(['(figure save skipped: ' name ')']);
end
end
