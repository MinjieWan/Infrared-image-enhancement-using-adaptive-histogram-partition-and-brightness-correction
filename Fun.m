function [ Fitness ] = Fun( D_thr,img,T,h_in,c_s,h_g,AHV,M0 )
% disp(['D_thr: ',num2str(D_thr)]);
count_s = size(AHV(:,1),1);
D = zeros(count_s,1);           
tmp = 0;
GCDF = zeros(count_s,1);
L = 256;
ga1 = 0.6;
for i = 1 : count_s
    if AHV(i,2) == 0 %background
        D(i,1) = (AHV(i,4) - AHV(i,3) + 1) * (L - D_thr ) / T;
        tmp = tmp + D(i,1);
    else
        t = [AHV(i,3),AHV(i,4)];
        ga = Gamma(t, h_in, ga1);
        GCDF(i,1) = c_s(i,1) ^ ga;
    end  
end

for i = 1 : count_s
    if AHV(i,2) == 1
        D(i,1) = (L  - D_thr - tmp) * GCDF(i,1) / sum(GCDF);
    end
end
D = round(D);

X  = zeros(count_s,2); 
for i = 1 : count_s
    if i == 1
        D0 = D_thr;
    end
    X(i,1) = D0;
    X(i,2) = D0 + D(i,1)-1 ;
    D0 = D0 + D(i,1);
end

y = zeros(256,1);
for i = 1 : count_s     
    if AHV(i,2) == 0
        I = (AHV(i,3):AHV(i,4))';
        y(AHV(i,3):AHV(i,4)) = D(i) * (I - AHV(i,3)) / (AHV(i,4) - AHV(i,3)) + X(i,1); 
    else
        hc = h_g(AHV(i,3):AHV(i,4));
        if sum(hc) == 0
            hc = zeros(AHV(i,4) - AHV(i,3) + 1, 1);
        else
            hc = hc / sum(hc);
        end
        cc = cumsum(hc);
        y(AHV(i,3):AHV(i,4)) = X(i,1) + (X(i,2) - X(i,1)) .* cc;
    end
       
end

F = double(img);
[row,col,~] =size(img);
for i = 1 : row
    for j = 1 : col
        F(i,j) = y(img(i,j)+1,1);
    end
end

Fitness = abs(mean2(F) - M0);



end

