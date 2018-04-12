% DongKyu Kim
% this function does pulse shaping and undoes pulse shaping
% I will work on this when time allows
% if trigger is 0 it's forward direction, if not then reverse.
function [out] = pulseshape(x,trigger)
	switch trigger
		case 0 % transmit
			out = x; % do nothing
		case 1 % receive
			out = x./abs(x); % normalize to let's say... account for noise
    end
end
% I think pulse shaping as a step closer to being realistic. However
% for this purpose, I think it's not necessary to do pulse shaping
% in the most advanced way, because what I'm going to do is shape the
% pulse and get it down to digital using the same method, which just
% adds noise to the signal. I considered doing this step, but this adds
% some bit delays to the message, and I don't think I want to modify
% my demodulator function which is already sort of problemaitc. But I'll
% leave this function alive so that I know what I need to do.