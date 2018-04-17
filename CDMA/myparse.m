function [ out ] = myparse(in, PN_SEQ,H)
% This is my parse function
% Function converts data into a human readable form ASCII characters
% Function XORs the input with the PN_sequence and despreads with the
% Hadamard transform.
frames = reshape(in,255,[]);
not_last = frames(:,1:end-1);
last_frame = frames(:,end); % in case last frame has fewer than 3 characters
not_last_PN = zeros(192,size(not_last,2));
for i = 1:size(not_last,2)
    not_last_PN(:,i) = xor(myBPSK(not_last(find(not_last(:,i)),i),'d'),PN_SEQ(1:192));
end
    temp = myBPSK(last_frame(find(last_frame)),'d');
    last_frame_PN = xor(temp,PN_SEQ(1:length(temp)));
    data = [not_last_PN(:);last_frame_PN(:)];
    decoded = myBPSK(H*reshape(myBPSK(data,'e'),8,[])/8,'d');
    ascii_dec = bi2de(reshape(decoded,8,[]).');
    out = [char(ascii_dec.')];
end

