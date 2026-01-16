function [U,S,V] = svdecon(X)
% Cheaper version of SVD


[m,n]       = size(X);
if  m <= n
    C       = X*X';
    [U,D]   = eig(C);
    clear C;

    [d,ix]  = sort(abs(diag(D)),'descend');
    U       = U(:,ix);

    if nargout > 2
        V   = X'*U;
        s   = sqrt(d);
        V   = bsxfun(@(x,c)x./c, V, s');
        S   = diag(s);
    end
else
    C       = X'*X;
    [V,D]   = eig(C);
    clear C;

    [d,ix]  = sort(abs(diag(D)),'descend');
    V       = V(:,ix);

    U       = X*V; 
    s       = sqrt(d);
    U       = bsxfun(@(x,c)x./c, U, s');
    S       = diag(s);
end
