function yHatPoint = lwlrPoint(point, xMat, yMat, k)
     [row , ~] = size(xMat);
     weights = zeros(row, row);
     for i = 1:1:row
         diffMat = point - xMat(i, :);
         weights(i,i) = exp(diffMat * diffMat.' / (-2.0 * (k ^ 2)));    
     end
     xTx = xMat.' * (weights * xMat);   
     if det(xTx) == 0
          disp('This matrix is singular, cannot do inverse');
     end
     theta = xTx^-1 * (xMat.' * (weights * yMat));          
     yHatPoint = point * theta;
end

