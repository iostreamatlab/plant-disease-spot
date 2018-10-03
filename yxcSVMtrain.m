function [alphaStar, bStar, SVIndex] = yxcSVMtrain(X, Y, C, kernel, sigma)


[num, dim] = size(X);
if dim ~= 2
    return;
end
if num ~= length(Y)
    return;
end
Y = Y(:);

H = (Y*Y').*yxcSVMkernel(X, X, kernel, sigma);
f = -ones(num, 1);
A = zeros(1, num);
b = 0;
Aeq = Y';
beq = 0;
lb = zeros(num, 1);
ub = C .* ones(num, 1);
x0 = zeros(num,1);
qp_options = optimset('MaxIter',10^3, 'LargeScale', 'off', 'Display','off');
[alphaStar,fval,exitflag,output] = quadprog(H,f,[],[],Aeq,beq,lb,ub,x0,qp_options);

nearZero = 10^-12;

alphaStar(find(abs(alphaStar) < nearZero)) = 0;
alphaStar(find(alphaStar > C - nearZero)) = C;

SVIndex = find(alphaStar > 0);

SVNotOnBoundIndex = find(alphaStar > 0 & alphaStar < max(alphaStar));
 

if ~isempty(SVNotOnBoundIndex)
    bStar = sum( Y(SVNotOnBoundIndex) -  H(SVNotOnBoundIndex, SVIndex) * alphaStar(SVIndex) .* Y(SVNotOnBoundIndex) )...
                 / length(SVNotOnBoundIndex);
else
    bStar = 0;
end
