function [ out ] = myBPSK(in, flag)
% This is my BPSK encoder / decoder
% Function converts in (bits) into BPSK signal or the opposite depending on
% the flag. 'e' if encoding, 'd' if decoding

switch flag
    case 'e'
        out = in*2-1;
    case 'd'
        out = in > 0;
end

end

