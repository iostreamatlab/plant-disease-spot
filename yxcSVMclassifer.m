function [YClassified, Z, err] = yxcSVMclassifer(Xtrain, Xtest, Y, alphaStar, bStar, kernel, sigma)



[num, dim] = size(Xtrain);
if dim ~= 2
    return;
end
if num ~= length(Y)
    return;
end
Y = Y(:);
alphaStar = alphaStar(:);

% Algorithm 5.4.12, step (4)
H = yxcSVMkernel(Xtrain, Xtest, kernel, sigma)';
Z = H * (alphaStar .* Y);
Z = Z + bStar;

YClassified(find(Z > 0)) = 1;
YClassified(find(Z <= 0)) = -1;
YClassified = YClassified';
if(size(Xtest) == size(Xtrain))
    if(Xtest == Xtrain)
        err = length(find(YClassified ~= Y)) / length(Y);
    end
else
    err = NaN;
end