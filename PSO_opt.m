function [ Thr ] = PSO_opt(img,T,h_in,c_s,h_g,AHV,M0, N, ger,limit)

d = 1;                         
vlimit = [-5, 5];              
w = 0.8;                       
c1 = 0.5;                    
c2 = 0.5;                     
x = zeros(N,d);
for i = 1:d  
    R = rand(N,1);
    x(:,i) = round(limit(i, 1) + (limit(i, 2) - limit(i, 1)) * R);
end 
v = rand(N, d);                 
xm = x;                         
ym = zeros(1, d);               
fxm = 100000*ones(N, 1);                 
fym = 100000;                      

% update 
iter = 1;  
% record = zeros(ger);             
% record_iter = zeros(ger,d);
fx = zeros(N,1);

while iter <= ger 
    
     for i = 1 : N
         fx(i,1) = Fun(x(i,:)', img,T,h_in,c_s,h_g,AHV,M0);
     end

     for i = 1:N    
        if fxm(i,1) > fx(i,1)  
            fxm(i,1) = fx(i,1);     
            xm(i,:) = x(i,:);       
        end   
     end  
     
    if fym > min(fxm)  
            [fym, nmax] = min(fxm);     
            ym = xm(nmax, :);         
    end
    
%     [record_iter(iter,1),nmax] = min(fx); 
%     record_iter(iter,2) = x(nmax, 1);

    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x); % velocity update
    
    for i = 1 : N
        if v(i,1) < vlimit(1,1)
            v(i,1) = vlimit(1,1);
        end
        if v(i,1) > vlimit(1,2)
            v(i,1) = vlimit(1,2);
        end
    end
    
    x =  x + v ; % position update  
    
    for i = 1 : N
        if x(i,1) < limit(1,1)
            x(i,1) = limit(1,1);
        end
        if x(i,1) > limit(1,2)
            x(i,1) = limit(1,2);
        end
    end
   
%     record(iter) = fym; 
    iter = iter+1;  
end  

Thr = round(ym);
end

