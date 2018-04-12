% DongKyu Kim
% this function undoes what Table 110 on p.45 does.
% or Table 108 p.44 column 1:2
function [out] = QPSKdecoding(x)
    x = mod(x,2*pi);
	out = zeros(length(x),2);
	for i = 1:length(x)
		if (x(i)>=pi/4) && (x(i)<3*pi/4)
			out(i,:) = [0,1];
		elseif (x(i) >= 3*pi/4) && (x(i)<5*pi/4)
			out(i,:) = [1,1];
		elseif (x(i) >= 5*pi/4) && (x(i)<7*pi/4)
			out(i,:) = [1,0];
		else
			out(i,:) = [0,0];
		end
	end
end	