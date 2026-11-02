%% circular_patch_theory  Analytical circular-patch radius (paper Eqs. 1-2).
%% Pure MATLAB/Octave. Reproduces the design equation behind the A = 14.5 mm
%% patch and the radius parametric study (paper Fig. 12) analytically; the
%% S11-vs-radius curves themselves need the full-wave solver.
close all

er = 4.4;                 % FR-4 relative permittivity
h_mm = 1.6;               % substrate thickness (mm)
h = h_mm / 10;            % -> cm (Eq. 1 uses cm)

% Effective permittivity, paper Eq. (2), with W = patch diameter = 2a.
ereff_of = @(a_cm) (er+1)/2 + (er-1)/2 * (1 + 10*h./(2*a_cm)).^0.5;

% Patch radius in cm from Eq. (1), iterated because ereff needs a.
radius_of = @(fr) local_radius(fr, er, h, ereff_of);

%% Design point: what radius resonates at 2.82 GHz?
fr0 = 2.82e9;
a0 = radius_of(fr0);
fprintf('Analytical radius at %.2f GHz: %.2f mm (paper: A = 14.5 mm)\n', ...
    fr0/1e9, a0*10);

%% Radius-vs-frequency curve (2.0 - 3.5 GHz).
fr = linspace(2.0e9, 3.5e9, 61);
a = zeros(size(fr));
for k = 1:length(fr)
    a(k) = radius_of(fr(k));
end
figure('Name', 'Circular patch design curve');
plot(fr/1e9, a*10, 'k-', 'LineWidth', 2);
hold on; grid on;
plot(fr0/1e9, a0*10, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot([2.0, 3.5], [14.5, 14.5], 'b--');
title('Circular patch radius vs resonant frequency (FR-4, h = 1.6 mm)');
xlabel('frequency / GHz'); ylabel('patch radius a / mm');
legend('paper Eqs. (1)-(2)', 'design point 2.82 GHz', 'paper choice A=14.5 mm', ...
    'Location', 'best');

%% Invert the curve: predicted resonance for the paper's radius sweep.
fprintf('\nRadius sweep (analytical fres vs paper Fig. 12 radii):\n');
fprintf('  a (mm) | predicted fres (GHz)\n');
for aTest = [13.5, 14.0, 14.5, 15.0, 15.5]
    fTest = interp1(a*10, fr, aTest, 'linear', 'extrap');
    fprintf('  %6.2f | %8.3f\n', aTest, fTest/1e9);
end

function a = local_radius(fr, er, h, ereff_of)
F = 8.791e9 / (fr*sqrt(er));          % cm
a = F;                                % seed
for it = 1:6                          % fixed-point on ereff(a)
    ereff = ereff_of(a);
    a = F / sqrt(1 + (2*h/(pi*ereff*F))*(log(pi*F/(2*h)) + 1.7726));
end
end
