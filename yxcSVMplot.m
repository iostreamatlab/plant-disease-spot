function yxcSVMplot(X, Y, SVIndex, alphaStar, bStar, kernel, sigma, plotAxis)


[num, dim] = size(X);
if dim ~= 2
    return;
end
if num ~= length(Y)
    return;
end
Y = Y(:);

plusIndex = find(Y==1);
minusIndex = find(Y==-1);

figure;
plot(X(plusIndex, 1), X(plusIndex, 2), 'r+');
if nargin == 8 
    axis(plotAxis);
end
hold on;
plot(X(minusIndex, 1), X(minusIndex, 2), 'g.');
hold on;
plot(X(SVIndex, 1), X(SVIndex, 2), 'ko');
hold on;

[cx, cy] = meshgrid(min( X(:, 1)) : .1 : max(X(:, 1)), min( X(:, 2)) : .1 : max(X(:, 2)));
[ncx, ncy] = size(cx);
Xtest = [cx(:) cy(:)];
[YClassified, Z]=yxcSVMclassifer(X, Xtest, Y, alphaStar, bStar, kernel, sigma);

Z = reshape(Z, ncx, ncy);


[C,h] = contour(cx, cy, Z, [0,0]);
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
hold on;
[C,h] = contour(cx, cy, Z, [1,1]);
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
hold on;
[C,h] = contour(cx, cy, Z, [-1,-1]);
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));