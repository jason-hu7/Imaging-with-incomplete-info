%% Problem set 1 question 2.2
clear, clc

n = 100;
H = shaw(n);
t = linspace(-pi/2, pi/2, n);

% In the first example, let's consider a linear s
s_1 = 2 * t;
y_1 = H * s_1';
figure()
plot(t, y_1, 'r'); hold on;
plot(t, s_1, 'b');
legend('y(t)', 's(t)');
xlabel('t'); ylabel('y(t) / s(t)');
title('Linear data')

% in the second example, let's consider an sinusoidal s
s_2 = sin(t) + cos(t);
y_2 = H * s_2';
figure()
plot(t, y_2, 'r'); hold on;
plot(t, s_2, 'b');
legend('y(t)', 's(t)');
xlabel('t'); ylabel('y(t) / s(t)');
title('sinusoidal data')

% in the second example, let's consider an exponential s
s_3 = exp(t);
y_3 = H * s_3';
figure()
plot(t, y_3, 'r'); hold on;
plot(t, s_3, 'b');
legend('y(t)', 's(t)');
xlabel('t'); ylabel('y(t) / s(t)');
title('exponential data')

%% question 2.3
% n = 20
H_20 = shaw(20);
spectrum_20 = svd(H_20);
normalized_s_20 = spectrum_20 / spectrum_20(1);

% n = 200
H_200 = shaw(200);
spectrum_200 = svd(H_200);
normalized_s_200 = spectrum_200 / spectrum_200(1);

% n = 2000
H_2000 = shaw(2000);
spectrum_2000 = svd(H_2000);
normalized_s_2000 = spectrum_2000 / spectrum_2000(1);

figure()
semilogy(normalized_s_20);
legend('N = 20');
xlabel('k'); ylabel('log of spectrum');
title('Spectrum');

figure()
semilogy(normalized_s_200);
legend('N = 200');
xlabel('k'); ylabel('log of spectrum');
title('Spectrum');

figure()
semilogy(normalized_s_2000);
legend('N = 2000');
xlabel('k'); ylabel('log of spectrum');
title('Spectrum');








