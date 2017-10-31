function [sest,VR] = GenLinInvSqrtQR(y,H,R,X,Q,VRform)
%  [sest,VR] = GenLinInvSqrtQR(y,H,R,X,Q,VRform)
%  
%       General Linear Estimation.
%       Calculates best estimates and covariance matrix of an
%       unknown vector using linear inverse problem
%       using qr decomposition
%       Intended for moderate-size problems
%
%  y - vector of measurements
%  H - observation matrix
%  R - Covariance matrix for measurement error
%         if scalar v, then R = v*eye(n) 
%  X - Drift matrix for s
%  Q - Covariance matrix for unknowns
%  sest  - Best estimates (posterior mean) for unknowns
%  VR:  if sixth argument is 'cov', VR'*VR = cov matrix
%       if sixth argument is 'inv', VR'*VR = inverse of cov matrix
%
%   PKK, June 2015, last tested July 2017
%

[m,p]=size(X);
[n,m1] = size(H);
[m2,m3] = size(Q);
[n1,n2] = size(R);

if var([m,m1,m2,m3]>0)
    error('check dimension m')
end

%Factorize Q
C = FactorQXMatrix2C(Q,X);
[KN,mw] = size(C);

%factorize R
if n1==n & n2==n 
    [B1,att1] = chol(R);
    if att1~=0
        error('R not positive definite')
    end
    B = inv(B1'); clear B1
    %test = norm(B'*B-inv(R))    
elseif n1==1 & n2==1
    B = 1/sqrt(R);
else
    error('R has wrong dimensions')
end

%solve LS problem
[sest, VR] = SolveLSQR([B*H;C],[B*y;zeros(KN,1)]);

tf = istriu(VR);
if tf ~=logical(1)
    error('VR not triangular')
end

if VRform == 'cov'
    opts.UT = true;
    VR = linsolve(VR,eye(m),opts);
    VR = VR';
elseif VRform == 'inv'
else
    error('error in sixth argument')
end

end
