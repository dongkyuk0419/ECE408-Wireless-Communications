%%
% DongKyu Kim
% Simple CDMA
% Wireless Communications
clc; clear all; close all;
%% Parameters and Data
Rcvd = load('Rcvd_Kim.mat');
Rcvd = Rcvd.Rcvd; % Received Signal
CR = 1e6; % Chip Rate
OS = 4; % Oversample
cpf = 255; % chip per frame
% RRCFR = 0.75; % RRC Filter Rolloff
B_RCOS = [0.0038,0.0052,-0.0044,-0.0121,-0.0023,0.0143,0.0044,-0.0385,-0.0563,0.0363,0.2554,0.4968,0.6025,0.4968,0.2554,0.0363,-0.0563,-0.0385,0.0044,0.0143,-0.0023,-0.0121,-0.0044,0.0052,0.0038].'; % Transmit Filter Coefficients
PN_SEQ = pngen([8,7,2,1],[1]); % PN Sequence % I flipped the polynomial because you did!
PN_SEQ_B = myBPSK(PN_SEQ,'e');
H = hadamard(8); % 8-ary Hadamard transform

% RRC Filter Receive
Rcvd_filt = filter(B_RCOS,1,Rcvd);
% Resolve Oversample
Rcvd_filt_down = downsample(Rcvd_filt,OS);

% Pilot Generation
pilot = zeros(8*4,1);
pilot_BPSK = myBPSK(pilot,'e');
pilot_hadamard = (pilot_BPSK*H(1,:)).';
pilot_rcvd = xor(myBPSK(pilot_hadamard(1:end-1),'d'),PN_SEQ(:).');
pilot_rcvd_fin = myBPSK(pilot_rcvd,'e');
% Well to be honest this could be one liner

% Figure out where the data starts and match pilot
frame_start = find(abs(Rcvd_filt_down(1:cpf))<0.1, 1, 'last' );
Pilot = Rcvd_filt_down(1+frame_start:cpf+frame_start);
Pilot2 = Rcvd_filt_down(1+frame_start+6*cpf:end);
% the magic number 6 came from
% length(Rcvd)/OS/cpf-1
hold on
figure(); plot(Pilot,'.'); % this looks reasonablly BPSK
k = myBPSK(real(Pilot),'d'); % there's an obvious boundary in the middle
% confirm = sum(pilot_rcvd == k) == cpf; % confirm that pilot matches.

figure; plot(Pilot2,'r.');
k2 = myBPSK(real(Pilot2),'d'); % there's an obvious boundary in the middle
k2 = not(k2);% actually for the last 255 the symbols are rotated 180 degrees
% confirm2 = sum(pilot_rcvd(1:length(k2)) == k2) == length(k2); % confirm that pilot matches.

% Extract phase rotation and frequency offset
k_phase = phase(Pilot);
k2_phase = phase(Pilot2);

freq_offset = (mean((k_phase(1:length(k2))-k2_phase)))/5;
% frequency offset is 0.3358 radians per 255 chips.
% So it would be 0.3358/2/pi/255*1e6 = 209.5Hz
freq_offset_per_chip = freq_offset / (cpf);

temp = k_phase - not(k)*pi;
temp2 = k2_phase - not(k2)*pi;

for i = 1: length(temp)-1
    store1(i) = temp(i+1) - temp(i);
    store2(i) = temp2(i+1) - temp2(i);
end
mean([store1,store2])

Data = Rcvd_filt_down(1+cpf+frame_start:frame_start+6*cpf);
Data_freq_shifted = Data.*exp(-1i*0.0028*[cpf:cpf+length(Data)-1]).*exp(1i*pi/4);
Rcvd_filt_down_freq_shifted = Rcvd_filt_down(1+frame_start:end).*exp(-1i*0.0028*[0:length(Rcvd_filt_down(1+frame_start:end))-1])*exp(1i*(pi/4+0.0028));
figure;plot(Rcvd_filt_down_freq_shifted(1:255),'.');
figure;
plot(Data_freq_shifted,'.');
axis([-2.5,2.5,-2.5,2.5])
figure;

plot(Rcvd_filt_down_freq_shifted,'.');
plot(Rcvd_filt_down(1:255),'r.'); hold on;plot(Rcvd_filt_down(256:256+255),'b.'); 

decoded = Data_freq_shifted - repmat(PN_SEQ_B.',1,5)
fin = myBPSK(real(decoded),'d');

