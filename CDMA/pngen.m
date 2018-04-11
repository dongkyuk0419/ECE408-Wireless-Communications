function [ seq ] = pngen( polynomial, initial_state)
% This is PN sequence generator
% polynomial is in the form suggested by the instruction
% ex) x4 + x3 + x2 +1 would be [4,3,2]
% initial_state is in the same form as the polynomial
order = polynomial(1);
state = zeros(order,1);
for i = 1:length(initial_state)
    state(initial_state(i)) = 1;
end
seq = zeros(2^order-1,1);

for i = 1:length(seq)
    seq(i) = state(end);
    temp = mod(sum(state(polynomial)),2);
    state(2:end) = state(1:end-1);
    state(1) = temp;
end

end

