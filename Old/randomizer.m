function M = randomizer(x,n,m)
M = zeros(m,n*length(x));
for j = 1:m
    temp = [];
    for i = 1:n
        temp = [temp x(randperm(length(x)))];
    end
    M(j,:) = temp;
end
end