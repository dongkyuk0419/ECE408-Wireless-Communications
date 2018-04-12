%DongKyu Kim
% clamper function. This function does what QPSKdecoding does (sort of) but
% instead of decoding it just clamps it to a phase value.
function [out] = clamper(x)
    [m,n] = size(x);
    out = x;
    for i = 1:m*n
            if (x(i)>=pi/4) && (x(i)<3*pi/4)
                out(i) = pi/2;
            elseif (x(i) >= 3*pi/4) && (x(i)<5*pi/4)
                out(i) = pi;
            elseif (x(i) >= 5*pi/4) && (x(i)<7*pi/4)
                out(i) = 3*pi/2;
            else
                out(i) = 0;
            end
    end
end
