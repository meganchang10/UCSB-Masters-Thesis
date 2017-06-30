function f = order(error_1,error_inf)

% This function tests the order of accuracy of a test. NOTE: The grid
% resolution must be doubling (i.e. 81x81, 161x161, 321x321) or this order
% is not valid

len = length(error_1)-1;
order_1 = zeros(len,1);
order_inf = zeros(len,1);

for i = 1:len
order_1(i) = log2(error_1(i+1)/error_1(i));
order_inf(i) = log2(error_inf(i+1)/error_inf(i));
end

order_1
order_inf

end