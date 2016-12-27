close all;
clear all;
%User inputs for defining range and velocity
R = input('Enter the distance of the target (5-50, units are km) ');
if R<5 || R>50
    error('That distance is not in range ');
end
R = R*1e3;
V = input('Enter the velocity of the target (m/s)');
if V<17.3 || V>170
    error('That velocity is not in range ');
end
%Determining Carrier frequency and PRF
c = 3e8; %speed of light
pulse_bw = c/(2*R); %pulse bandwidth
pulse_width = 1/pulse_bw ;%pulse width
npi = 10; %number fof pulses to integrate (lowers SNR)
f = 13e9;
wc = 2*pi*f;
bound1  = c/(2*R);
bound2 = (abs(V)*4*f)/c;
mur = 50e3;%maximum unambigous range
% %iur = 5e3;%minimum unambigous range
PRF = (abs(bound1-bound2)/2)+max([bound1 bound2]);
if PRF > mur
    PRF = mur;
end
%PRF= 8e3%test value
%CARRIER SIGNAL
t = 0:1/(20*f):2e-9;
fc = sin(wc*t);
a = min([length(fc) length(t)]);
% figure(1)
% plot(t(1:a),fc(1:a));
% title('Carrier Signal');
% xlabel('time');
% ylabel('amplitude');

%PULSES
ts = 0:1/f:(1/f)*(5e6);
A = 1/PRF;
on = (pulse_width)/A;
off = 1 - on;
Tr = [on*A,off*A];
pulse = double(mod(ts,sum(Tr))<Tr(1));
%  figure(2);
%  plot(ts,pulse);
%  title('Pulses');
%  xlabel('time'); 
%  ylabel('amplitude');

%re-time scale carrier signal and pulses
 fc = interp(fc,64);
 fc = [fc fc fc fc fc fc fc];
 t = 0:1/(10*f):7e-5;
 pulse = downsample(pulse,5);
%  figure(3)
%  a = min([length(t) length(pulse)]);
%  plot(t(1:a),pulse(1:a))
%  hold on;
%  a = min([length(t) length(fc)]);
%  plot(t(1:a),fc(1:a));
%  title('Scaled and Truncated');
%  xlabel('time'); 
%  ylabel('amplitude');

%GENERATE TRANSMITTED SIGNAL
 a = min([length(t) length(pulse) length(fc)]);
 fo = fc(1:a).*pulse(1:a);
%  figure(4);
 a = min([length(t) length(fo)]);
%  plot(t(1:a),fo(1:a))
%  title('Transmitted signal');
%  xlabel('time'); 
%  ylabel('amplitude');
  
% Doppler shift
 lamda = c/f;
 fd = (-2*V)/lamda;
 wd = 2*pi*fd;
 tp = 0:abs(1/(3*(f+fd))):(on)*A; 
 shift = sin((wc+wd)*tp); 
%  figure(5);
%  plot(tp,shift);
%  title('Doppler Shifted Pulse');
%  xlabel('time'); 
%  ylabel('amplitude');

%Add non-pulsed portion (basically a bunch of zeros)
 add =floor((length(shift)/on)*off);
 shift = [shift, zeros(1,add)];
 step = tp(2)-tp(1);
 t= 0:step:length(shift)*step;
 a = min([length(t)  length(shift)]);
%  figure(6);
%  plot(t(1:a),shift(1:a));
%  title('Doppler Shifted Signal');
%  xlabel('time'); 
%  ylabel('amplitude');

% time delay for received signal
td = (2*R)/c;
add = floor(td/step);
shift = [zeros(1,add) shift];
%increase t
t = 0:step:length(shift)*step;
% figure(7);
a = min([length(t)  length(shift)]);
% plot(t(1:a),shift(1:a));
% title('Received Signal w/ Time delay');
% xlabel('time'); 
% ylabel('amplitude');

threshold = ones(1,length(shift)); %used for determining threshold value

%ADD GAUSSIAN NOISE
SNR =input('Enter desired signal to noise ratio');
noise = awgn(shift,SNR,'measured','linear');
threshold = awgn(threshold,SNR,'measured','linear');
threshold = max(threshold);
figure(8);
a = min([length(t) length(noise)]);
plot(t(1:a),noise(1:a))
title('Gaussian White Noise + Received Signal');
xlabel('time'); 
ylabel('amplitude');

%estimation of range based on Time delay of received signal
i = find(noise>0.5);
td = t(i(1));
est_Range = (c/2)*td
signal= shift(i(1):i(end));
% figure(13);
% plot(abs(fft(signal)))%for debugging
% title('fft plot of received signal');



% %correlation
y = xcorr(shift,fc);
a = min([length(t),length(y)]);
figure(10)
plot(t(1:a),y(1:a));
title('Correlation');
xlabel('time'); 
ylabel('amplitude');

%autocorrolation
y = xcorr(signal);
a = min([length(t),length(y)]);
figure(11)
plot(t(1:a),y(1:a));
title('Autocorrelation');
xlabel('time'); 
ylabel('amplitude');
% %Attempt at matched filtering
shift = fft(signal);
shift = downsample(shift,2);
fo = fft(fo);
shift = shift(1:length(fo));
matched = shift.*fo;
figure(12);
plot(abs(matched));
title('Frequency Domain Matched Filter Response');
xlabel('frequency'); 
ylabel('amplitude');
matched = ifft(matched);
a = min([length(t) length(matched)]);
figure(9);
plot(t(1:a),matched(1:a))
title('Matched Filter response');
xlabel('time'); 
ylabel('amplitude');
 
