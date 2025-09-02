clear all; clc;

N = 16; % number of antennas
theta0 = 0;
% ============ load "signal.mat" file =============
load('signal.mat')
% ====================== end ======================

% ======= compute steering vector of theta0 =======

% ====================== end ======================

Pn = db2pow(-10);
Ns = 10000;

% compute the MVDR beamformer
wmvdr = mvdr(x, theta0);
% compute the MSINR beamformer
a_theta0 = exp(1j * pi * (0:N-1)' * sind(theta0));
Thetas = -90 : 0.1 : 89.9;
a = steervec(N, Thetas);  % 先建好 a
wmsinr = msinr(a_theta0, x, sp, a, Thetas);

% compute the MMSE beamformer
wmmse = mmse(x, sp);

% plot
Thetas = -90 : 0.1 : 89.9;
a = steervec(N, Thetas);
bpmvdr = wmvdr' * a;
bpmsinr = wmsinr' * a;
bpmmse = wmmse' * a;

%normalize
% Compute power
P_mvdr = abs(bpmvdr).^2;
P_msinr = abs(bpmsinr).^2;
P_mmse = abs(bpmmse).^2;

% Normalize to max = 1
P_mvdr = P_mvdr / max(P_mvdr);
P_msinr = P_msinr / max(P_msinr);
P_mmse = P_mmse / max(P_mmse);


figure;
semilogy(Thetas, P_mvdr, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MVDR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)

figure;
semilogy(Thetas, P_msinr, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MSINR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)

figure;
semilogy(Thetas, P_mmse, 'linewidth', 1.5)
axis([-90, 90, 1e-8, 1e2])
title('MMSE', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('P(\theta)', 'fontsize', 11)


figure;
plot(Thetas, angle(bpmvdr))
title('MVDR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)
figure;
plot(Thetas, angle(bpmsinr))
title('MSINR', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)
figure;
plot(Thetas, angle(bpmmse))
title('MMSE', 'fontsize', 14)
xlabel('angle', 'fontsize', 11)
ylabel('\angle P(\theta)', 'fontsize', 11)