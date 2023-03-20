function out = num2logidx_ycgosu(in)
% useful with [~, ~, in] = unique(numidx);

tempidx = zeros(size(in, 1), max(in));

for i = 1:max(in)
    tempidx(:, i) = (in == i);
end

out = logical(tempidx);

end