function J_lim = estimate_limit(p, p_range)

J_lim = 0;

n = length(p);
JJ = zeros(1,n);
for i = 1:n
    if p(i) >= p_range(1,i) && p(i) <= p_range(2,i)
        JJ(i) = 1;
    end
end
if sum(JJ) < n
    return
end

J_lim = 1;


end