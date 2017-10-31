clear, clc
%% Canonical straight-line tomography problem
% consider 2-D image with 8 rows and 12 columns(96 pixels)

%build the y vector(observation)
y = zeros(20,1);
y(1) = 12;
y(9:20) = 1;

%build the canonical straightline tomography data(first row straight line)
s = zeros(96,1);
for i = 1:8:96
    s(i) = 1;
end

%let's set up the H observation matrix
H = zeros(20, 96);
for i = 1:1:8
    for j =i:8:96
        H(i, j) = 1;
    end
end

k = 1;
for i = 9:1:20
    for j = k:1:k+7
        H(i,j) = 1;
    end
    k = k + 8;
end

H_rank = rank(H);
y_hat = H * s;
s_hat = H' * inv(H * H') * y;

%% 1.2
spectrum_h = svd(H);
[U, S, V] = svd(H);
normalized_spectrum_h = spectrum_h / spectrum_h(1);
figure()
semilogy(normalized_spectrum_h);
xlabel('k'); ylabel('log of spectrum');
title('Spectrum');

%% 1.3
figure()
for i = 1:1:20
    subplot(4,5,i);plot(H(i,:));
end

figure()
for i = 1:1:20
    subplot(4,5,i);plot(V(:,i));
end

%% 1.5
figure()
imagesc(H*H');
colorbar
figure()
imagesc(H'*H);
colorbar
