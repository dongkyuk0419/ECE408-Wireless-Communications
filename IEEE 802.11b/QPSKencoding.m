% DongKyu Kim
% this function does what Table 110 on p.45 does.
% or Table 108 p.44 column 1:2
function [out] = QPSKencoding(x)
    out=zeros(length(x(:,1)),1);
    for i=1:length(out)
        switch char(x(i,:))
            case char([0,0])
                out(i)=0;
            case char([0,1])
                out(i)=pi/2;
            case char([1,1])
                out(i)=pi;
            case char([1,0])
                out(i)=-pi/2;
        end
    end
end

