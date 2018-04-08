% This Matlab file demomstrates an IR image enhancement algorithm based on Minjie Wan et al's paper:
% "Infrared image enhancement using adaptive histogram partition and brightness correction" published by Remote Sensing (2008)
% Author: Minjie Wan, all rights reserved
% E-mail: minjiewan1992@njust.edu.cn 
% http://www4.comp.polyu.edu.hk/~cskhzhang/

clear all;
close all;
clc;

% parameter setting
sig = 0.7;                        
window = 3;
sigma = 3;                      
thr_hat = 0;
L = 256;
ga1 = 0.6;  
alpha = 0.005;
W = 11;
N = 10; % particle number of PSO
ger = 10; % iteration of PSO

% load the input IR image
path = '1.bmp';
img = imread(path);
[row,~,dim] = size(img);
if dim ~= 1
    img = rgb2gray(img);
end
img = img(1:row-3,:);
figure(1);
imshow(img);
title('input IR image');

tic;
[row,col] = size(img);
h_in = imhist(img);
h_in(256,1) = 0;
h_in = h_in / sum(h_in);

C = zeros(row,col); % calculate the local contrast weighted distribution
for i = 2 : row - 1
    for j = 2 : col - 1
        m1 = abs(img(i,j) - img(i+1,j)); 
        m2 = abs(img(i,j) - img(i-1,j));
        m3 = abs(img(i,j) - img(i,j+1));
        m4 = abs(img(i,j) - img(i,j-1));
        C(i,j) = double((m1+m2+m3+m4)/4);
    end
end

h_g = zeros(256,1);
for i = 2 : row - 1
    for j = 2 : col - 1 
        tmp = img(i,j) + 1;
        h_g(tmp,1) = h_g(tmp,1) + Local_contrast([i,j], C);
    end
end
h_g = h_g / sum(h_g);

%LOWESS
count = (1:256)';
y = h_in;                          
x = [ones(length(y),1) count];
y_hat = lwlr(x, x, y, sigma);      
y_hat = y_hat';
y_hat (y_hat< thr_hat) = 0;

%local minima examination
flag = 0;
count = 1;
while flag == 0                  % the first non-zero
    if h_in(count) == 0
        count = count + 1;
    else
        flag = 1;
        m1 = count; %
    end
end

flag = 0;
count = 256;
while flag == 0                    % the last non-zero
    if h_in(count) == 0
        count = count - 1;
    else
        flag = 1;
        m2 = count;         
    end
end

% s:stors all the valleys
s = zeros(1,1);
s(1,1) = m1;
s(1,2) = h_in(m1);
count = 2;
for i = m1 + (W + 1)/2  : m2 - (W + 1)/2 
    if y_hat(i) == min(y_hat(i - (W - 1)/2  : i + (W - 1)/2))
        s(count,1) = i;
        s(count,2) = h_in(i);
        count = count + 1;
    end
end
s(count,1) = m2;
s(count,2) = h_in(m2);

% fore-and background recognition
count_s = count - 1;                % number of intervals
c_s = zeros(count_s,1);
AHV = zeros(count_s,4);
for i  = 1 : count_s
    if i == 1
        AHV(i,3) = s(i,1); 
    else 
        AHV(i,3) = s(i,1) + 1;
    end
    AHV(i,4) = s(i+1,1); 
    t = AHV(i,4) - AHV(i,3) + 1;   
    sub_h = h_in(AHV(i,3):AHV(i,4));
    c_s(i,1) = sum(sub_h);
    AHV(i,1) = c_s(i,1) * (row*col) / t; 
end
AHV(:,1) = AHV(:,1) / max(AHV(:,1));
level = graythresh(AHV(:,1));     % Otsu

for i = 1 : count_s
    if AHV(i,1) >= level          
        AHV(i,2) = 0;             % background 
    else
        AHV(i,2) = 1;             % foreground
    end
end
s1 = 0;
for i = 1 : count_s
    if AHV(i,2) == 0
        s1 = s1 + AHV(i,4) - AHV(i,3);
    end
end

T = m2 - m1 + 1;                  
path = 'Ref.tif';
Ref = imread(path);
[~,~,dim] = size(Ref);
if dim ~= 1
    Ref = rgb2gray(Ref);
end
M0 = mean2(Ref);

% brightness correction using PSO
s1 = 0;
for i = 1 : count_s
    if AHV(i,2) == 0
        s1 = s1 + AHV(i,4) - AHV(i,3);
    end
end
limit = [0, 255-s1];  % position limitation of particles
D_thr  = PSO_opt( img,T,h_in,c_s,h_g,AHV,M0, N, ger,limit );

% calculate the enhanced image using the optimized minimal grayscale D_thr
D = zeros(count_s,1);            
tmp = 0;
GCDF = zeros(count_s,1);
for i = 1 : count_s
    if AHV(i,2) == 0
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
for i = 1 : row
    for j = 1 : col
        F(i,j) = y(img(i,j)+1,1);
    end
end
eps = abs(mean2(F)-mean2(Ref));

% display the enhanced image
F = F / max(max(F));
figure(2);
imshow(F);
title('enhanced image');
toc
