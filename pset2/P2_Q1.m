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
s_hat = pinv(H) * y;
[U, S, V] = svd(H);
Vr = V(:, 1:rank(H));
Vn = V(:, (rank(H)+1):size(V,2));

%% after the matrix is built. We will dive into peoble 1 of problemset2
% case 1
% in the first case we have actual image given by value 1 at every pixel
% except at pixels on the second row where it is 2.
s_1 = ones(96, 1);
for i = 2:8:96
    s_1(i) = 2;
end
% figure()
% imagesc(reshape(s_1, [8,12]));colorbar;

y_1 = [12;24;12;12;12;12;12;12;9;9;9;9;9;9;9;9;9;9;9;9];
s_hat_1 = pinv(H) * y_1;

figure()
subplot(1,2,1)
imagesc(reshape(s_hat_1, [8,12]));colorbar;title('estimate');
subplot(1,2,2)
imagesc(reshape(s_1, [8,12]));colorbar;title('true');

err1 = immse(s_1, s_hat_1);
fprintf('\n The mean-squared error is %0.4f\n', err1);

%compute the rsquare for s_1
R_square_1 = (norm(Vr*Vr'*s_1)^2)/(norm(s_1)^2);

%% case 2
% in the second case the image is 10 at every pixel except at column 3
% where the value is 0
s_2 = ones(96, 1);
s_2(1:96) = 10;
s_2(17:24) = 0;
% figure()
% imagesc(reshape(s_2, [8,12]));colorbar;
y_2 = [110;110;110;110;110;110;110;110;80;80;0;80;80;80;80;80;80;80;80;80];
s_hat_2 = pinv(H) * y_2;

figure()
subplot(1,2,1)
imagesc(reshape(s_hat_2, [8,12]));colorbar;title('estimate');
subplot(1,2,2)
imagesc(reshape(s_2, [8,12]));colorbar;title('true');

err2 = immse(s_2, s_hat_2);
fprintf('\n The mean-squared error is %0.4f\n', err2);

% compute the r square value for s_2
R_square_2 = (norm(Vr*Vr'*s_2)^2)/(norm(s_2)^2);

%% case 3
% in the thrid case the image is 1 everywhere except at pixels
% 1,9,10,18,19,27,28,36,37,45,46,54,55,63,64,72 where it is 0
s_3 = ones(96, 1);
for i = 1:9:64
    s_3(i) = 0;
end
for i = 9:9:72
    s_3(i) = 0;
end
% figure()
% imagesc(reshape(s_3, [8,12]));colorbar;

y_3 = zeros(20,1);
y_3(1:8) = 10;
y_3(9) = 7;
y_3(10:16) = 6;
y_3(17) = 7;
y_3(18:20) = 8;
s_hat_3 = pinv(H) * y_3;

figure()
subplot(1,2,1)
imagesc(reshape(s_hat_3, [8,12]));colorbar;title('estimate');
subplot(1,2,2)
imagesc(reshape(s_3, [8,12]));colorbar;title('true');

err3 = immse(s_3, s_hat_3);
fprintf('\n The mean-squared error is %0.4f\n', err3);

% compute the r square value for s_3
R_square_3 = (norm(Vr*Vr'*s_3)^2)/(norm(s_3)^2);

%% case 4
% the image is 1 everywhere except at pixels 1,2,3,9,10,11,17,18,19 where
% it is 1
s_4 = ones(96,1);
s_4(1:3) = 0;
s_4(9:11) = 0;
s_4(17:19) = 0;
% figure()
% imagesc(reshape(s_4, [8,12]));colorbar;

y_4 = zeros(20,1);
y_4(1:3) = 9;
y_4(4:8) = 12;
y_4(9:11) = 5;
y_4(12:20) = 8;
s_hat_4 = pinv(H) * y_4;

figure()
subplot(1,2,1)
imagesc(reshape(s_hat_4, [8,12]));colorbar;title('estimate');
subplot(1,2,2)
imagesc(reshape(s_4, [8,12]));colorbar;title('true');

err4 = immse(s_4, s_hat_4);
fprintf('\n The mean-squared error is %0.4f\n', err4);

% compute the r square value for s_2
R_square_4 = (norm(Vr*Vr'*s_4)^2)/(norm(s_4)^2);

%% p1.3
% Take columns 20 through 24 of V and multiply each of these 5 columns by a
% random variable from 0 to 1
rng(0,'twister');
a=0;b=1;
r = (b-a).*rand(5,1) + a;
sum_v20_24 = zeros(96,1);
for i = 20:24
    sum_v20_24 = sum_v20_24 + r(i-19)*V(:,i);
end

figure()
imagesc(reshape(sum_v20_24, [8,12]));colorbar;