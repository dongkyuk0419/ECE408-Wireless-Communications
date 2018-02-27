% DongKyu Kim
% ECE408 - Wireless Communication
% Professor Hoerning
% sequence generator
% This function generates BPSK modulated signal of size N.
% Currently only supports BPSK (M = 2)

function [out] = seqgen(M,N)
M = 2; % force BPSK
x = randi([0,M-1],N,1); % generate bit sequence
out = x*2-1; %BPSK modulation
end