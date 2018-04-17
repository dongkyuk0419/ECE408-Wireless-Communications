function [ out ] = myBPSK(in, flag)
% This is my BPSK encoder / decoder
% Function converts in (bits) into BPSK signal or the opposite depending on
% the flag. 'e' if encoding, 'd' if decoding, 'ds' for special decoding for
% to simulate switch

switch flag
    case 'e'
        out = in*2-1;
    case 'd'
        out = in > 0;
    case 'ds'
        out = zeros(size(in));
        for i = 1:length(in)
            if in(i) > 0.5
                out(i) = 1;
            elseif abs(in(i))< 0.5
                out(i) = 0;
            elseif in(i) < -0.5
                out(i) = -1;
            end
        end
end

end

