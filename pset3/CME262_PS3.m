clear, clc
% CME 262 Problem set 3, we need to compute and plot the posterior mean and
% the 95% credibility interval in each of the sceanario
load HW3P1.txt
y = HW3P1;  % y is a 400 x 1 vector
m = 400;
H = shaw(m);

%% Q1. constant mean and exponential covariance
% l = pi/2 and variance = 100, variance of observation is 1e-8
% to use GenlinInv, get y, H, R, X, Q
%  y - vector of measurements ***
%  H - observation matrix ***
%  R - Covariance matrix for measurement residual
%  X - Drift matrix for s, must have rank >0
%  Q - Covariance matrix for unknowns
l = pi/2;
var = 100;

R0 = 1e-8*eye(m);
expQ = @(h) var*exp(-h/l);

%construct covariance matrix for unkowns
Q0 = zeros(m,m);
t = linspace(-pi/2,pi/2,m);
for i = 1:m
    for j = 1:m
        Q0(i,j) = expQ(norm(t(i) - t(j)));
    end
end

X0 = ones(m,1); %drift matrix for s
[s_0,V_0,LAMBDA_0,MU_0] = GenLinInv(y,H,R0,X0,Q0);

[s_1,VR] = GenLinInvSqrtQR(y,H,R0,X0,Q0, 'cov');

figure()
plot(t, s_0, 'r', t, s_1, '--', 'linewidth', 2);
legend('s using GenLinInv', 's using GenLinInvSqrtQR');
xlabel('t');
ylabel('s_0 / s_1');

%% Q2 replacing the exponential with exponential squared aka Gaussian
GaussQ = @(h) var*exp(-(h^2)/(l^2));
Q1 = zeros(m,m);
for i = 1:m
    for j = 1:m
        Q1(i,j) = GaussQ(norm(t(i) - t(j)));
    end
end
[s_2,V_2,LAMBDA_2,MU_2] = GenLinInv(y,H,R0,X0,Q1);

[s_3,VR3] = GenLinInvSqrtQR(y,H,R0,X0,Q1, 'cov');

figure()
plot(t, s_2, 'r', t, s_3, '--', 'linewidth', 2);
legend('s using GenLinInv', 's using GenLinInvSqrtQR');
xlabel('t');
ylabel('s_2 / s_3');

%% Q3 try one more model Nugget

Q2 = var .* eye(m);

[s_4,V_4,LAMBDA_4,MU_4] = GenLinInv(y,H,R0,X0,Q2);

[s_5,VR5] = GenLinInvSqrtQR(y,H,R0,X0,Q2, 'cov');

figure()
plot(t, s_4, 'r', t, s_5, '--', 'linewidth', 2);
legend('s using GenLinInv', 's using GenLinInvSqrtQR');
xlabel('t');
ylabel('s_4 / s_5');

%% Q5 generate 5 conditional realizations for q1 and q2

s_u = simv(Q0, X0);
v = zeros(5,m);
s_c = zeros(5,m);
for i = 1:5
    v(i,:) = mvnrnd(zeros(m,1), R0);
    s_c(i,:) = s_u + LAMBDA_0 * (y + v(i,:)' - H*s_u);
end

figure()
plot(s_c(1,:), 'linewidth', 1, 'DisplayName', '1st conditional realization'); hold on;
plot(s_c(2,:), 'linewidth', 1, 'DisplayName', '2nd conditional realization'); hold on;
plot(s_c(3,:), 'linewidth', 1, 'DisplayName', '3rd conditional realization'); hold on;
plot(s_c(4,:), 'linewidth', 1, 'DisplayName', '4th conditional realization'); hold on;
plot(s_c(5,:), 'linewidth', 1, 'DisplayName', '5th conditional realization'); hold on;
legend('show')
