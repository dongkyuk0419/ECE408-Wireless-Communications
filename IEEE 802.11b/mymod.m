%DongKyu Kim
% modulates signal based on the rate
% 1,2Mbits/s is done by barker sequence spreading
% 5.5, 11Mbits/s is done by CCK (complementary code keying) modulation
function [out] =  mymod(x,rate)
barker = [1 -1 1 1 -1 1 1 1 -1 -1 -1].'; % 11 chip barker sequence (p.42)
C = @(p1,p2,p3,p4) [exp(1i*(p1+p2+p3+p4)),exp(1i*(p1+p3+p4)), exp(1i*(p1+p2+p4)),-exp(1i*(p1+p4)),exp(1i*(p1+p2+p3)),exp(1i*(p1+p3)),-exp(1i*(p1+p2)),exp(1i*(p1))];
% CCK code word (p.43)

switch rate
    case 1 %Barker code spreading with DBPSK (1 bit / symbol) (p.42)
        modulator = comm.DBPSKModulator;
        s = step(modulator,x);
        out = reshape(barker*s.',[],1); %transmit left most barker chip first then go on.
    case 2 %Barker code spreading with DQPSK (2 bits / symbol) (p.42)
        modulator = comm.DQPSKModulator(0,'BitInput',true, 'SymbolMapping', 'Gray');
        % 0 phase shift,  assume bit inputs, and gray coding (p.43)
        if mod(length(x),2)~=0
            x = [x;0];
        end
        s = step(modulator,x);
        out = reshape(barker*s.',[],1);      
    case 5.5 %CCK (4 bits/symbol) (p.44)
        if mod(length(x),4)~=0
            x = [x;zeros(4 - mod(length(x),4),1)];
        end
        x = reshape(x,4,[]).';
        p1 = QPSKencoding(x(:,1:2));
        p1(1:2:end) = p1(1:2:end) + pi; % odd symbols offset by pi
        % (p.45)
        p2 = x(:,3)*pi+pi/2;
        p3 = zeros(size(p1));
        p4 = x(:,4)*pi;
        s = C(p1,p2,p3,p4).';
        out = s(:); %align symbol order correctly
        
    case 11 %CCK (8 bits/symbol) (p.45)
        if mod(length(x),8)~=0
            x = [x;zeros(8 - mod(length(x),8),1)];
        end        
        x = reshape(x,8,[]).';
        p1 = QPSKencoding(x(:,1:2));
        p1(1:2:end) = p1(1:2:end) + pi; % odd symbols offset by pi
        % Table 110 (p.45)
        p2 = QPSKencoding(x(:,3:4));
        p3 = QPSKencoding(x(:,5:6));
        p4 = QPSKencoding(x(:,7:8));
        s = C(p1,p2,p3,p4);
        out = reshape(s.',[],1); %align symbol order correctly

    otherwise
        disp('invalid rate');
end