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
PN_SEQ = pngen([8,7,2,1],1); % PN Sequence % I flipped the polynomial because you did!
H = hadamard(8); % 8-ary Hadamard transform

% RRC Filter Receive
Rcvd_filt = filter(B_RCOS,1,Rcvd)   ;
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
num_frames_till_last = length(Rcvd)/OS/cpf-1;
frame_start = find(abs(Rcvd_filt_down(1:cpf))<0.1, 1, 'last' );
Pilot = Rcvd_filt_down(1+frame_start:cpf+frame_start);
Pilot2 = Rcvd_filt_down(1+frame_start+num_frames_till_last*cpf:end);
% figure(); plot(Pilot,'.'); % this looks reasonablly BPSK
% plot(Pilot2,'r.');
k = myBPSK(real(Pilot),'d'); % there's an obvious boundary in the middle
% confirm = sum(pilot_rcvd == k) == cpf; % confirm that pilot matches.

k2 = myBPSK(real(Pilot2),'d'); % there's an obvious boundary in the middle
k2 = not(k2);% actually for the last 255 the symbols are rotated 180 degrees
% confirm2 = sum(pilot_rcvd(1:length(k2)) == k2) == length(k2); % confirm that pilot matches.

% Extract phase rotation and frequency offset
k_phase = phase(Pilot);
k2_phase = phase(Pilot2);

% freq_offset = (mean((k_phase(1:length(k2))-k2_phase)))/5;
% freq_offset_per_chip = freq_offset / (cpf);
% Actually this is not correct. I don't really get why.

temp = k_phase - not(k)*pi;
temp2 = k2_phase - not(k2)*pi;

store1 = zeros(1,length(temp)-1);
store2 = zeros(1,length(temp2)-1);
for i = 1: length(temp)-1
    store1(i) = temp(i+1) - temp(i);
end
for i = 1: length(temp2)-1
    store2(i) = temp2(i+1) - temp2(i);
end

freq_offset_per_chip = median([mod(store1,2*pi),mod(store2,2*pi)]);
% frequency offset is 0.0028 radians per chip.
% So it would be 1e6 * 0.0028 / 2 / pi, 445Hz
phase_offset = median(angle(Pilot.*exp(-1i*freq_offset_per_chip*(0:length(Pilot)-1))) - not(k)*pi);

% Data Frequency unshift
Data = Rcvd_filt_down(1+cpf+frame_start:frame_start+num_frames_till_last*cpf);
Data_freq_shifted = Data.*exp(-1i*freq_offset_per_chip*(cpf:cpf+length(Data)-1)).*exp(-1i*phase_offset);

% Despread and decode
Data_no_pilot = Data_freq_shifted - repmat(pilot_rcvd_fin,1,num_frames_till_last-1);
decoded1 = myBPSK(real(Data_no_pilot),'ds');
string = myparse(decoded1,PN_SEQ,H(6,:))
