%% Assignment 1.1
clear all
clc
numerator = [0.505 1.01 0.505];
denominator = [1 0.7478 0.2722];
% Linear time-invariant filter (Butterworth filter)
F = tf(numerator,denominator,1,'Variable','q^-1');
plotoptions = bodeoptions;
plotoptions.FreqScale = 'linear';
bode(F,plotoptions)
%% Assignment 1.2
t = -50:1:50;
% Input r is a linear ramp.
r = t;
[u,y] = assignment_sys_18(r);
plot(t,u)
title('Input u(t) over time with signal r(t)=t')
xlabel('tiempo [s]')
ylabel('u')
%% Assignment 2.1
Range = [-3 3];
Band = [0 0.7003];
N = 1024;
SINEDATA = [128,128,1];
% Input r is a sum of sines at 128 frequencies
r = idinput(N,'SINE',Band,Range,SINEDATA);
% plot power spectrum with periodic signal input
rfft = fft(r);
Fs = 0.7*60; % Sampling frequency in rad/hour
L = length(r);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(rfft(1:L/2))) % Plot until Nyquist frequency
title('Single-Sided Amplitude Spectrum of r(t)')
xlabel('Frequency (rad/hour)')
ylabel('Amplitude')
%pwelch(r)
%cpsd(r,r);
%pspectrum(r,'FrequencyResolution',0.07)
%plot(1:1:1024,r)
% Sample cross covariance R.^(N)yu(tao)
%RyueN = cra(data,128);
%title('Sample cross covariance R.^Nyu(tao)') %Choose Hamming window lag

%% Assignment 2.2
figure;
[u,y] = assignment_sys_18(r);
data = iddata(y,u,1);
%Gsmooth = etfe(data,100,128); % Hamming window lag = 20
Gnon = etfe(data,[],128); % no Hamming window lag
%bode(Gnon,Gsmooth)
bode(Gnon);
%legend('No smoothing','\gamma=100')
legend('No smoothing')
%% Assignment 2.3
[G,wout] = freqresp(Gnon);
v = fft(y) - G.*fft(u);

%vfft = fft(Rv);
Fs = 0.007; % Sampling frequency in rad/hour
L = length(v);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(vfft(1:L/2))) % Plot until Nyquist frequency
title('Magnitud of Spectrum of v(t)=y(t)-G(q,\theta_N)*u(t)')
xlabel('Frequency (rad/hour)')
ylabel('Amplitude')
%% Assignment 3.1
%PRBS 
%% Assignment 3.2
Range = [-3,3];
Band = [0 0.7];
r = idinput(3000,'prbs',Band,Range);
% plot power spectrum with periodic signal input
rfft = fft(r);
Fs = 0.7*60; % Sampling frequency in rad/hour
L = length(r);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(rfft(1:L/2))) % Plot until Nyquist frequency
title('Single-Sided Amplitude Spectrum of r(t)')
xlabel('Frequency (rad/hour)')
ylabel('Amplitude')

%% Assignment 4.1
[u,y] = assignment_sys_18(r);


