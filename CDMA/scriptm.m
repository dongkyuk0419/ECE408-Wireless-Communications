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
RRCFR = 0.75; % RRC Filter Rolloff
B_RCOS = [0.0038,0.0052,-0.0044,-0.0121,-0.0023,0.0143,0.0044,-0.0385,-0.0563,0.0363,0.2554,0.4968,0.6025,0.4968,0.2554,0.0363,-0.0563,-0.0385,0.0044,0.0143,-0.0023,-0.0121,-0.0044,0.0052,0.0038].'; % Transmit Filter Coefficients
PN_SEQ = pngen([8,7,2,1],[1]); % PN Sequence % I flipped the polynomial because you did!
H = hadamard(8); % 8-ary Hadamard transform

