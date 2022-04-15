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
Band = [0 0.7];
N = 1024;
SINEDATA = [128,128,1];
% Input r is a sum of sines at 128 frequencies
[r freq] = idinput(N,'SINE',Band,Range,SINEDATA);
freq(1);
% plot power spectrum with periodic signal input
rfft = fft(r);
Fs = 2*pi; 
L = length(r);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(rfft(1:L/2))) % Plot until Nyquist frequency
title('Single-Sided Amplitude Spectrum of r(t)')
xlabel('Frequency (rad/s)')
ylabel('Amplitude')

%% Assignment 2.2
figure;
[u,y] = assignment_sys_18(r);
data = iddata(y,u,1);
% Frequency response estimation
Gnon = etfe(data,[],128); % no Hamming window lag
%Gsmooth = etfe(data,30,128); % Hamming window lag = 60
bode(Gnon);
legend('No smoothing')
%% Assignment 2.3
% Estimated spectrum
[spec_yu,F] = cpsd(y,u,[],0,128);
[spec_yy,F] = cpsd(y,y,[],0,128);
[spec_uu,F]= cpsd(u,u,[],0,128);
spec_vv = spec_yy - (abs(spec_yu).^2)./spec_uu;
% Plot estimated spectrum
semilogx(F,20*log(spec_vv),'-*')
hold on
semilogx(F,20*log(spec_yy),'-*')
hold off
title('Estimated noise spectrum')
ylabel('Magnitude (dB)')
xlabel('Frequency (rad/s)')

%% Assignment 3.1
%PRBS 
%% Assignment 3.2
Range = [-3,3];
Band = [0 0.5];  
r = idinput(3000,'prbs',Band,Range);
[u,y] = assignment_sys_18(r);
% plot power spectrum with periodic signal input
rfft = fft(u);
Fs = 2*pi; % Sampling frequency in rad/s
L = length(u);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(rfft(1:L/2))) % Plot until Nyquist frequency
title('Single-Sided Amplitude Spectrum of u(t)')
xlabel('Frequency (rad/s)')
ylabel('Amplitude')
% Bode plot of the model G
figure;
z0 = iddata(y,u);
G0 = etfe(z0,[],3000); % no Hamming window lag
bode(G0,'.',{0.01,pi})
%% Assignment 4.1
NN = struc(3,4,1);
M = oe(z0,NN);
% Bode plot of system G and the model G0(OE)
figure;
bode(G0,'y.',M,'b',{0.1,pi});
legend('No smoothing',strcat('OE','[',num2str(NN),']'))
figure;
resid(z0,M,'corr');
[ys,fit,ic] = compare(z0,M)
ysim = idsim(u,M);
RMSE = sqrt(sum((y(:)-ysim(:)).^2) / numel(y) ); % Root Mean Squared Error
MAE = sum(abs(ysim(:)-y(:))) / numel(y); % Mean Absolute Error
%% Assignment 5.1
Range = [-3,3];
Band = [0 0.7]; 
r = idinput(1000,'rbs',Band,Range);
[u,y] = assignment_sys_18(r);
% plot power spectrum with periodic signal input r(t)
rfft = fft(u);
Fs = 2*pi; % Sampling frequency in rad/s
L = length(u);
w = (0:L-1)*Fs/L;
stem(w(1:L/2),abs(rfft(1:L/2))) % Plot until Nyquist frequency
title('Single-Sided Amplitude Spectrum of u(t)')
xlabel('Frequency (rad/s)')
ylabel('Amplitude')
% Monte Carlo simulation:
number_of_runs = 100;
B_m = [];
F_m = [];
cov_m = [];
for n=1:number_of_runs
    [u ,y] = assignment_sys_18(r);
    data_m = iddata(y,u,1);
    sys = oe(data_m,NN);
    [A,B,C,D,F,dA,dB,dC,dD,dF] = polydata(sys);
    cov_data = getcov(sys);
    cov_m = [cov_m; cov_data];
    sys_g = append(sys);
    B_m = [B_m; B];
    F_m = [F_m; F];
end
%% Assignment 5.3
th_var = cov_data
mean(B_m);
mean(F_m);
Bsim_stat = var(B_m)
Fsim_stat = var(F_m)
%% Assignment 5.4
B_init = median(B_m)
F_init = median(F_m)
M_init = idpoly([],B_init,[],[],F_init);
Bq = [];
Fq = [];
for n=1:number_of_runs
    [u ,y] = assignment_sys_18(r);
    data_m = iddata(y,u,1);
    M_oe = oe(data_m,M_init);
    [A,B_q,C,D,F_q,dA,dB,dC,dD,dF] = polydata(M_oe);
    covq = getcov(M_oe);
    Bq = [Bq; B_q];
    Fq = [Fq; F_q];
end
th_var = cov_data
mean(Bq);
mean(Fq);
Bsim_stat = var(Bq)
Fsim_stat = var(Fq)
%% Assignment 6.1
figure;
z0 = iddata(y,u);
G0 = etfe(z0,[],1000); % no Hamming window lag
bode(G0,'.',{0.01,pi})
%%
NN_bj = struc(4,2,3,4,1);
M_bj = bj(z0,NN_bj);
% Bode plot of system G and the model G0(OE)
figure;
bode(G0,'y.',M_bj,'b',{0.1,pi});
legend('No smoothing',strcat('BJ','[',num2str(NN_bj),']'))
figure;
%opt = residOptions('MaxLag',20);
resid(z0,M_bj,'corr');
%[y,fit,ic] = compare(data,M_bj)
th_var = getcov(M_bj);
th_var(:,5:9) = [];
th_var(5:9,:) = []
