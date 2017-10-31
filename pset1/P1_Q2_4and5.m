clear,clc
%% q2.4
% for m = n = 400, plot the most important vectors that span the rowspace
% of H.
H_400 = shaw(400);
rank_H_400 = rank(H_400);
[U,S,V] = svd(H_400);
figure()
for i = 1:1:20
    subplot(4,5,i);plot(V(:,i));
end

%% q2.5
n = 400;
H = shaw(n);
t = linspace(-pi/2, pi/2, n)';
s = 2 * exp(-6 *(t - 0.8).^ 2) + exp(-2 * (t + 0.5).^ 2);
y = H * s;
s_hat = pinv(H)*y;
figure()
plot(t, s, 'r', t, s_hat, '--', 'linewidth', 2);
legend('real s', 'estimated s');
xlabel('t');
ylabel('s / s_hat');

%% q2.6
v1 = normrnd(0,1e-12,n,1);
y_noise1 = H * s + v1;
s_noise1 = pinv(H) * y_noise1;
figure();
subplot(3,1,1);
plot(t, s, 'r', t, s_noise1, '--', 'linewidth', 2);
legend('real s', 'estimated s with normal noise 1e-12');
xlabel('t');
ylabel('s / s_hat');

v2 = normrnd(0,1e-15,n,1);
y_noise2 = H * s + v2;
s_noise2 = pinv(H) * y_noise2;
subplot(3,1,2);
plot(t, s, 'r', t, s_noise2, '--', 'linewidth', 2);
legend('real s', 'estimated s with normal noise 1e-15');
xlabel('t');
ylabel('s / s_hat');

v3 = normrnd(0,1e-8,n,1);
y_noise3 = H * s + v3;
s_noise3 = pinv(H) * y_noise3;
subplot(3,1,3);
plot(t, s, 'r', t, s_noise3, '--', 'linewidth', 2);
legend('real s', 'estimated s with normal noise 1e-8');
xlabel('t');
ylabel('s / s_hat');

