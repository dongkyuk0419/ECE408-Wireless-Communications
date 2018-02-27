% DongKyu Kim
% ECE408 - Wireless Communication
% Professor Hoerning
% Alamouti Coding
% Script

clc; clear all; close all;

% Parameters
M = 2; %BPSK modulation flag
k = 100; % number of monte carlo runs
N = 1024; % message length

for i = 1:k
    bit_stream = seqgen(M,N);
    
end
