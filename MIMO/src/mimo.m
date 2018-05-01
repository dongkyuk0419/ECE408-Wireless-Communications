%Vishnu Kaimal
%DongKyu Kim
%MIMO/OFDM
%MIMO part

clc; clear all; close all;
%% params

numIter = 3;  % The number of iterations of the simulation
nSym = 5e1;    % The number of symbols per packet
M = 16;        % Modulation order
Mt = 2;     %number of transmitters
Mr = 2;     %number of receivers

EbNo = -10:2:30; %EbNo range to iterate over for plot
SNR_Vec = EbNo + 10*log10(log2(M)); %SNR conversion from EbNo
lenSNR = length(SNR_Vec); 
ber = zeros(3,numIter,lenSNR); %ber store
H = sqrt(1/2)*(randn(Mr,Mt,nSym*log2(M),numIter)+1j*randn(Mr,Mt,nSym*log2(M),numIter));%channel for each iteration


%% precoding

U = zeros(Mr,Mt,nSym*log2(M));
V = zeros(Mr,Mt,nSym*log2(M));
tx_chan = zeros(Mr,1,nSym*log2(M));
tx_process = zeros(Mr,1,nSym*log2(M));

for i = 1:numIter
    
    bits = randi([0 M-1],Mt, 1, nSym*log2(M));     % Generate random bits
    tx = qammod(bits,M);  % modulate the signal
        
    %construction of flat fading channel
    
    for k=1:nSym*log2(M)
        [U(:,:,k),~,V(:,:,k)] = svd(H(:,:,k,i));
        tx_chan(:,:,k) = H(:,:,k,i)*V(:,:,k)*tx(:,:,k);
    end
        
    for j = 1:lenSNR
         noise = sqrt(1/2)*(randn(Mr,1,nSym*log2(M))+1j*randn(Mr,1,nSym*log2(M)));
        txNoisy = tx_chan + 10^(-1*SNR_Vec(j)/20)*noise;
        for k=1:nSym*log2(M)
            tx_process(:,:,k) = U(:,:,k)'*txNoisy(:,:,k);
        end
        rx = qamdemod(tx_process,M);

        % Compute and store the BER for this iteration
        [~, ber(1,i,j)] = biterr(bits, rx); 

   end  
end

%% zero forcing


for i = 1:numIter
    
    bits = randi([0 M-1],Mt, 1, nSym*log2(M));     % Generate random bits
    tx = qammod(bits,M);  % modulate the signal
    
    %construction of flat fading channel
    
    for k=1:nSym*log2(M)
        tx_chan(:,:,k) = H(:,:,k,i)*tx(:,:,k);
    end
    
    
    for j = 1:lenSNR
 
        noise = sqrt(1/2)*(randn(Mr,1,nSym*log2(M))+1j*randn(Mr,1,nSym*log2(M)));
        txNoisy = tx_chan + 10^(-1*SNR_Vec(j)/20)*noise;

        for k=1:nSym*log2(M)
            W(:,:,k) = (H(:,:,k,i)'*H(:,:,k,i))^-1*H(:,:,k,i)';
            tx_process(:,:,k) = W(:,:,k)*txNoisy(:,:,k);
        end
        
        
        rx = qamdemod(tx_process,M);

        % Compute and store the BER for this iteration
        [~, ber(2,i,j)] = biterr(bits, rx); 

   end  
end

%% MMSE


for i = 1:numIter
    
    bits = randi([0 M-1],Mt, 1, nSym*log2(M));     % Generate random bits
    tx = qammod(bits,M);  % modulate the signal
    
    %construction of flat fading channel
    
    for k=1:nSym*log2(M)
        tx_chan(:,:,k) = H(:,:,k,i)*tx(:,:,k);
    end
    
    for j = 1:lenSNR
 
        noise = sqrt(1/2)*(randn(Mr,1,nSym*log2(M))+1j*randn(Mr,1,nSym*log2(M)));
        txNoisy = tx_chan + 10^(-1*SNR_Vec(j)/20)*noise;

        for k=1:nSym*log2(M)
            W(:,:,k) = (H(:,:,k,i)'*H(:,:,k,i)+eye(Mt)*10^(-1*SNR_Vec(j)/20))^-1*H(:,:,k,i)';
            tx_process(:,:,k) = W(:,:,k)*txNoisy(:,:,k);
        end
        
        rx = qamdemod(tx_process,M);

        % Compute and store the BER for this iteration
        [~, ber(3,i,j)] = biterr(bits, rx); 

   end  
end

%% plot

ber = mean(ber,2); %take mean across all iterations

%plotting ber for different equalization techniques
semilogy(EbNo,ber(1,:),'DisplayName','precoding');
hold on;
semilogy(EbNo,ber(2,:),'DisplayName','Zero-Forcing');
hold on;
semilogy(EbNo,ber(3,:),'DisplayName','MMSE');

title('BER plots for different MIMO Schemes');
xlabel('EbNo(dB)');
ylabel('Bit Error Rate');
legend('show');