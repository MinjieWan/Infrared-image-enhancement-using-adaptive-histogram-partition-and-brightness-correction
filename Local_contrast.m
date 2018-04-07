function [ fc ] = Local_contrast( pos, C )

i = pos(1,1);
j = pos(1,2);
c = C(i,j);
Jmax = max(max(C));

eps1 = 0.01;
eps2 = 0.99;

B = atanh(2 * eps1 - 1);
A = (atanh(2 * eps2 - 1) - B) / Jmax;

fc = 0.5 * (tanh(A * c + B) + 1);

end

