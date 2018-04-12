%DongKyu Kim
% demodulates signal based on the rate
% 1,2Mbits/s is done by barker sequence spreading
% 5.5, 11 Mbits/s is done by CCK (complementary code keying) demodulation 
% (I couldn't find the RAKE decoder that the standard refers to so I made my own)
function [out] =  mydemod(x,rate)
barker = [1 -1 1 1 -1 1 1 1 -1 -1 -1].'; % 11 chip barker sequence (p.42)
C = @(p1,p2,p3,p4) [exp(1i*(p1+p2+p3+p4)),exp(1i*(p1+p3+p4)), exp(1i*(p1+p2+p4)),-exp(1i*(p1+p4)),exp(1i*(p1+p2+p3)),exp(1i*(p1+p3)),-exp(1i*(p1+p2)),exp(1i*(p1))];
% CCK code word (p.43)

switch rate
    case 1 %Barker code spreading with DBPSK (1 bit / symbol) (p.42)
        unspread = reshape(x,11,[]);
        s = barker.'*unspread/11;
        demodulator = comm.DBPSKDemodulator;
        out = step(demodulator,s.');
    case 2 %Barker code spreading with DQPSK (2 bits / symbol) (p.42)
        unspread = reshape(x,11,[]);
        s = barker.'*unspread/11;
        demodulator = comm.DQPSKDemodulator(0,'BitOutput',true, 'SymbolMapping', 'Gray');
        out = step(demodulator,s.');
    case 5.5 %CCK (4 bits/symbol) (p.44)
		unspread = reshape(x,8,[]).';
		unwrapped = mod(angle(unspread),2*pi);
        unwrapped = clamper(unwrapped);
		% reverse what was done
        p1 = mod(unwrapped(:,8),2*pi);
        p2 = mod(unwrapped(:,5)-p1,2*pi);
        p4 = mod(unwrapped(:,2)-p1,2*pi);
% 		p1 = mod((unwrapped(:,6) + unwrapped(:,8))/2,2*pi);
% 		p4 = mod((unwrapped(:,2) + unwrapped(:,4)-pi)/2-p1,2*pi);
% 		p2 = mod(((unwrapped(:,1) + unwrapped(:,3)-pi)/2-p1-p4+(unwrapped(:,5)+unwrapped(:,7)-pi)/2-p1)/2,2*pi);
% Trying something like differential mode in electronics doesn't work...
        p1(1:2:end) = p1(1:2:end)-pi;
		out = zeros(4,length(p1));
		out(1:2,:) = QPSKdecoding(p1).';
		out(3,:) = (p2>pi); % decision becomes this
		out(4,:) = (mod(p4+pi/2,2*pi)>pi); %shift the decision boundary so that it's easier
		out = out(:);
	case 11 %CCK (8 bits/symbol) (p.45)
		unspread = reshape(x,8,[]).';
		unwrapped = mod(angle(unspread),2*pi);
		p1 = mod((unwrapped(:,8)),2*pi);
        p2 = mod((unwrapped(:,7)-pi-p1),2*pi);
        p3 = mod((unwrapped(:,6)-p1),2*pi);
        p4 = mod((unwrapped(:,4)-pi-p1),2*pi);
        p1(1:2:end) = p1(1:2:end)-pi;
        out = zeros(8,length(p1));
        out(1:2,:) = QPSKdecoding(p1).';
        out(3:4,:) = QPSKdecoding(p2).';
        out(5:6,:) = QPSKdecoding(p3).';
        out(7:8,:) = QPSKdecoding(p4).';
        out = out(:);

    otherwise
        disp('invalid rate');
end


