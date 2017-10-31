function y = simv(Q,X)

if nargin == 1
    m = length(Q);
    V = (chol(Q))';
    y = V*randn(m,1);
    
elseif nargin == 2
    [m,p] = size(X); [q,r] = qr(X);
    q1 = q(1:m, 1:p); q2 = q(1:m, p+1:m);
    Qa = q2'*Q*q2; V = (chol(Qa));
    ya = V * randn(m-p,1); y = q2 * ya;
end