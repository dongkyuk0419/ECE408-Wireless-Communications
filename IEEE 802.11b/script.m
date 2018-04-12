%% IEEE 802.11b
% PHY Layer
% DongKyu Kim
% ECE408 Wireless Communications

clear all; close all; clc;

%% Parameters
EbNo = -4:1:12;
n_packet = 100;    % number of packets sent
packet_size = 1e5; % number of bits sent
datarates = [1,2,5.5,11]; % all in Mbits/s (p.14)
bps = [1,2,4,8]; % bits per symbol | see mymod.m for references
%% Simulation
ber_store = zeros(length(datarates),length(EbNo));
for i = 1:length(datarates) % data rate loop
    for ii = 1:length(EbNo) % EbNo loop
        n_bits_tx = 0;
        n_bits_incor = 0;
        snr = EbNo(ii) + 10*log10(bps(i));
        
        for k = 1:n_packet % packet loop
            n_bits_tx = n_bits_tx + packet_size;
            m = randi([0,1],packet_size,1);
            m_tx = mymod(m,datarates(i));
			m_tx_p = pulseshape(m_tx,0);
            channel_tx = awgn(m_tx_p,snr,'measured');
			m_rx_p = pulseshape(channel_tx,1);
            m_rx = mydemod(m_rx_p,datarates(i));
            n_bits_incor = n_bits_incor + sum(m~=m_rx);
        end
        ber_store(i,ii) = n_bits_incor/n_bits_tx;
    end
end

%% Plotting

figure();
semilogy(EbNo, ber_store(1,:),'-o')
title('BER vs EbNo for different data rates in IEEE 802.11b');
xlabel('EbNo (dB)');
ylabel('Bit Error Rate');
hold on;

for i=2:4
    semilogy(EbNo,ber_store(i,:),'-o')
end
legend('1Mbps','2Mbps','5.5Mbps', '11Mbps');



