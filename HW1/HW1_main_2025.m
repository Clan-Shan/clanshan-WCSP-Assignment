clear all;
close all;

%% Parameters of Rayleigh model

BW = 1000;
Ts = 1/BW;
N_FFT = 1024;

%% Rayleigh model

ch_Rayleigh = HW1_Rayleigh_2025(Ts);

% Draw the power spectrum density for Rayleigh model
ch_Rayleigh_auto = autocorr(abs(ch_Rayleigh), 'NumLags', length(ch_Rayleigh)-1);
ch_Rayleigh_auto = [ch_Rayleigh_auto zeros(1,N_FFT-length(ch_Rayleigh))];
ch_Rayleigh_psd = fftshift(abs(fft(ch_Rayleigh_auto)));
f_idx = -BW/2:(BW/1024):(BW-1)/2;
H1 = figure(1);
plot(f_idx, ch_Rayleigh_psd);
ylabel('PSD');
xlabel('freq. (Hz)');
grid;

% Evaluation the PDF of the magitude fading coefficient for Rayleigh model
mag_ch_Rayleigh_realization = abs(ch_Rayleigh);
x = 0:0.01:8;
sigma2 = 0.5;
pdf_theory = (x/sigma2).*exp(-x.^2/(2*sigma2));
pdf_env = ksdensity(mag_ch_Rayleigh_realization,x);
H2 = figure(2);
plot(x,pdf_env,'-',x,pdf_theory,'*');
legend('Simulated','Theoretic');
xlabel('mag.');
ylabel('PDF of mag.');
grid;

% Evaluation the PDF of the phase fading coefficient for Rayleigh model
phase_ch_Rayleigh_realization = angle(ch_Rayleigh);
x = -pi:0.01:pi;
pdf_phase = ksdensity(phase_ch_Rayleigh_realization,x);
H3 = figure(3);
plot(x,pdf_phase,'-');
legend('Simulated');
xlabel('phase');
ylabel('PDF of phase');
axis([-pi pi 0 0.5]);
grid;

%% Parameters of Jakes model
nu = 12e4; %7e4
c = 3e8;
fc= 7e9; %3.5e9
f_max = (nu/3600)/((c)/(fc));
BW = 1000;
Ts = 1/BW;
M = 60;
N_FFT = 1024;
fprintf("%f\n", f_max)
%% Jakes model

ch_Jakes = HW1_Jakes_2025(M, f_max, Ts);

% Draw the power spectrum density for Jakes model
ch_Jakes_auto = autocorr(real(ch_Jakes), "NumLags", length(ch_Jakes)-1);
% ch_Jakes_auto = xcorr(ch_Jakes,'normalized');

% ch_Jakes_auto = ch_Jakes_auto/ch_Jakes_auto(length(ch_Jakes));
% ch_Jakes_auto = ch_Jakes_auto(length(ch_Jakes):2*length(ch_Jakes)-1);

ch_Jakes_auto = [ch_Jakes_auto zeros(1,N_FFT-length(ch_Jakes))];
ch_Jakes_psd = fftshift(abs(fft(ch_Jakes_auto)));



f_idx = -BW/2:(BW/1024):(BW-1)/2;
H4 = figure(4);
plot(f_idx, ch_Jakes_psd);
ylabel('PSD');
xlabel('freq. (Hz)');
grid;

ch_Jakes_psd = ch_Jakes_psd/trapz(f_idx, ch_Jakes_psd);
% plot(f_idx, ch_Jakes_psd);
% f_bar = trapz(f_idx, f_idx.*ch_Jakes_psd)/trapz(f_idx, ch_Jakes_psd);
f_bar = sum(f_idx.*ch_Jakes_psd)/sum(ch_Jakes_psd);
% f_RMS = sqrt(trapz(f_idx, ((f_idx-f_bar).^2).*ch_Jakes_psd)/trapz(f_idx, ch_Jakes_psd))
f_RMS = sqrt(sum(((f_idx-f_bar).^2).*ch_Jakes_psd)/sum(ch_Jakes_psd));
fprintf("f_RMS : %f\n", f_RMS)

Sd = @(f) (1/pi/f_max).*(1./sqrt(1-(f/f_max).^2));
Sd2 = @(f) f*(1/pi/f_max).*(1./sqrt(1-(f/f_max).^2));
Sd3 = @(f) ((f-f_bar).^2).*(1/pi/f_max).*(1./sqrt(1-(f/f_max).^2));
f_bar = integral(Sd2,-f_max,f_max)/integral(Sd,-f_max,f_max);
f_RMS2 = sqrt(integral(Sd3,-f_max,f_max)/integral(Sd,-f_max,f_max));
plot(f_idx, Sd(f_idx));
fprintf("f_RMS2 : %f\n", f_RMS2)

% Evaluation the PDF of the magitude fading coefficient for Jakes model
mag_ch_Jakes_realization = abs(ch_Jakes);
x = 0:0.01:6;
sigma2 = 0.5;
pdf_theory = (x/sigma2).*exp(-x.^2/(2*sigma2));
pdf_env = ksdensity(mag_ch_Jakes_realization,x);
H5 = figure(5);
plot(x,pdf_env,'-',x,pdf_theory,'*');
legend('Simulated','Theoretic');
xlabel('mag.');
ylabel('PDF of mag.');
grid;

% Evaluation the PDF of the phase fading coefficient for Jakes model
phase_ch_Jakes_realization = angle(ch_Jakes);
x = -pi:0.01:pi;
pdf_phase = ksdensity(phase_ch_Jakes_realization,x);
H6 = figure(6);
plot(x,pdf_phase,'-');
legend('Simulated');
xlabel('phase');
ylabel('PDF of phase');
axis([-pi pi 0 0.5]);
grid;
