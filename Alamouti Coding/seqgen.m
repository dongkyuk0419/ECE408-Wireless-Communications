% DongKyu Kim
% ECE408 - Wireless Communication
% Professor Hoerning
% sequence generator
% This function generates BPSK modulated signal of size N.
% Currently only supports BPSK (M = 2)

function [out1,out2] = seqgen(M,N)
M = 2; % force BPSK
out2 = randi([0,M-1],N,1); % generate bit sequence
out1 = out2*2-1; %BPSK modulation
end