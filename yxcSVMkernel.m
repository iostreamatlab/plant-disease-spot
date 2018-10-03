function K = yxcSVMkernel(a, b, kernel, sigma)

[n1, dim1] = size(a);
[n2, dim2] = size(b);

if dim1 ~= dim2
    error('columns of a and b must agree');
end

switch kernel
    case 'linear'
        K = a*b';
    case 'poly'
        K = (a * b' + 1).^sigma;
    case 'rbf'
        [N1, d] = size(a);
        [N2, d] = size(b);
        dist2 = repmat(sum((a.^2)', 1), [N2 1])' + ...
                repmat(sum((b.^2)',1), [N1 1]) - 2*a*(b');
        K = exp(-dist2/(sigma^2));        
    otherwise
        error('kernel type must be specified.');
end

