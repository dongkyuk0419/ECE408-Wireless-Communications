% DongKyu Kim
% ECE408 - Wireless Communication
% Professor Hoerning
% Alamouti Coding
% Script

clc; clear all; close all;

% Parameters
M = 2; % Symbol size
B = log2(M); % bits per symbol
N = 1e7; % message length
EbNo = 0:5:50;
snr = EbNo + 10*log10(B);

% Rayleigh Channel (MATLAB)
ts = 1e-5;
fd = 130;
rc_1 = rayleighchan(ts,fd);
rc_2 = rayleighchan(ts,fd);
rc_3 = rayleighchan(ts,fd);
rc_4 = rayleighchan(ts,fd);

% Data
[original, modulated] = seqgen(M,N);

% No Diversity
tx_temp = filter(rc_1,modulated);
tx = zeros(N,length(snr));
rx = zeros(N,length(snr));
for i = 1:length(snr)
    tx(:,i) = awgn(tx_temp,snr(i),'measured');
    rx(:,i) = real(tx(:,i))>0;
end
[~,ber_nd] = biterr(rx,original);




