function [ ga ] = Gamma( T,his, ga1 )

t1 = T(1,1);
t2 = T(1,2);
t = (t1:t2)';
h = his(t1:t2,1);
h = h / sum(h);
M = sum(t .* h);
ga = M / 255 * ga1;

end

