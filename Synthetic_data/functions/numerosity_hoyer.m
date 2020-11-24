function [num] = numerosity_hoyer(x)
% numerosity - Finds the numerosity of input matrix x
%        Each column of x is taken to be a vector and
%        any element greater than threshold is assumed
%        to be nonzero.
% hoyer sparse inspire sparsity when it reach 1;
%

[r, c] = size(x);
num = zeros(1, c);

for i = 1:c
    D=norm(x(:,i),1)/norm(x(:,i),2);
    D(isnan(D)) = 1;
    sparsehoyer=(sqrt(r)-D)/(sqrt(r)-1);
    sparsehoyer(isnan(sparsehoyer)) = 1;
    num(1, i) = sparsehoyer;%     
end
