clear,clc
%P2_Q2
%% first we need to formulate the covariance matrix
sigma_square = 0.8;
l = 1.5;

Q_ij = @(d_ij, sigma_square, l) sigma_square * exp(-d_ij./l);
d_ij = ones(96, 96);

for i = 1:96
    for j = 1:96
        if i <= 8
            if j <= 8
                d_ij(i,j) = sqrt((j-i)^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-48-i)^2 + 6^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-56-i)^2 + 7^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-64-i)^2 + 8^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-72-i)^2 + 9^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-80-i)^2 + 10^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-88-i)^2 + 11^2);
            end
        elseif (i>8 && i<=16)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-8-i)^2 + 2^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-16-i)^2 + 3^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-24-i)^2 + 4^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-32-i)^2 + 5^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-40-i)^2 + 6^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-48-i)^2 + 7^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-56-i)^2 + 8^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-64-i)^2 + 9^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-72-i)^2 + 10^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-80-i)^2 + 11^2);
            end
        elseif (i>16 && i<=24)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-48-i)^2 + 6^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-56-i)^2 + 7^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-64-i)^2 + 8^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-72-i)^2 + 9^2);
            end
        elseif (i>24 && i<=32)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-48-i)^2 + 6^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-56-i)^2 + 7^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-64-i)^2 + 8^2);
            end
        elseif (i>32 && i<=40)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-48-i)^2 + 6^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-56-i)^2 + 7^2);
            end
        elseif (i>40 && i<=48)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-48-i)^2 + 6^2);
            end
        elseif (i>48 && i<=56)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-40-i)^2 + 5^2);
            end
        elseif (i>56 && i<=64)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+56)^2 + 7^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-32-i)^2 + 4^2);
            end
        elseif (i>64 && i<=72)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+64)^2 + 8^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+56)^2 + 7^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-24-i)^2 + 3^2);
            end
        elseif (i>72 && i<=80)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+72)^2 + 9^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+64)^2 + 8^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+56)^2 + 7^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-8-i)^2 + 1^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-16-i)^2 + 2^2);
            end
        elseif (i>80 && i<=88)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+80)^2 + 10^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+72)^2 + 9^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+64)^2 + 8^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+56)^2 + 7^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            end
        elseif (i>88 && i<=96)
            if j <= 8
                d_ij(i,j) = sqrt((j-i+88)^2 + 11^2);
            elseif (j>8 && j<=16)
                d_ij(i,j) = sqrt((j-i+80)^2 + 10^2);
            elseif (j>16 && j<=24)
                d_ij(i,j) = sqrt((j-i+72)^2 + 9^2);
            elseif (j>24 && j<=32)
                d_ij(i,j) = sqrt((j-i+64)^2 + 8^2);
            elseif (j>32 && j<=40)
                d_ij(i,j) = sqrt((j-i+56)^2 + 7^2);
            elseif (j>40 && j<=48)
                d_ij(i,j) = sqrt((j-i+48)^2 + 6^2);
            elseif (j>48 && j<=56)
                d_ij(i,j) = sqrt((j-i+40)^2 + 5^2);
            elseif (j>56 && j<=64)
                d_ij(i,j) = sqrt((j-i+32)^2 + 4^2);
            elseif (j>64 && j<=72)
                d_ij(i,j) = sqrt((j-i+24)^2 + 3^2);
            elseif (j>72 && j<=80)
                d_ij(i,j) = sqrt((j-i+16)^2 + 2^2);
            elseif (j>80 && j<=88)
                d_ij(i,j) = sqrt((j-i+8)^2 + 1^2);
            elseif (j>88 && j<=96)
                d_ij(i,j) = sqrt((j-i)^2 + 0^2);
            end
        end
    end
end

%construct covariance matrix for unknowns using meshgrid
% [X,Y] = meshgrid(1:12, 1:8);
% [Xi, Xj] = meshgrid(X(:),X(:));
% [Yi, Yj] = meshgrid(Y(:), Y(:));
% d = sqrt((Xi-Xj).^2+(yi-yj).^2);
% Q_mesh = exp(-d/l);

Q = Q_ij(d_ij,sigma_square, l);
figure()
imagesc(Q); colorbar;

%% second let's construct the H matrix

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

% H_rank = rank(H);
% y_hat = H * s;
% s_hat = pinv(H) * y;
% [U, S, V] = svd(H);
% Vr = V(:, 1:rank(H));
% Vn = V(:, (rank(H)+1):size(V,2));

%% 2.2 
R = (10^-8) * eye(20);
prior_mean = 0;
% equatino 105
posterior_mean = prior_mean + inv((inv(Q) + H'*inv(R)*H))*H'*inv(R)*(y-H*prior_mean);
% equation 111
Xi = (H*Q*H'+R)\(y-H*prior_mean);
posterior_mean2 = prior_mean + Q * H' * Xi;

%% 2.3
Q_2 = Q_ij(d_ij, 1, 2);
R_2 = (10^-8) * eye(20);
prior_mean2 = 0;
lambda = ((H*Q_2*H'+R_2)\(H*Q_2))';
spectrum_lambda = svd(lambda);
[U_lambda, S_lambda, V_lambda] = svd(lambda);
normalized_spectrum_lambda = spectrum_lambda / spectrum_lambda(1);

figure()
semilogy(normalized_spectrum_lambda);
xlabel('k'); ylabel('log of spectrum');
title('Spectrum');

figure()
for i = 1:1:19
    subplot(4,5,i);plot(U_lambda(:,i));
end

%% 2.4
% real image
s_4 = ones(96,1);
s_4(1:3) = 0;
s_4(9:11) = 0;
s_4(17:19) = 0;
% data
y_4 = zeros(20,1);
y_4(1:3) = 9;
y_4(4:8) = 12;
y_4(9:11) = 5;
y_4(12:20) = 8;
y_4_noise = y_4 + normrnd(0,1e-8,20,1);

posterior_mean4 = lambda*(y_4_noise);
posterior_Q = -lambda*H*Q + Q;
figure()
imagesc(reshape(posterior_mean4, [8,12]));colorbar;

%% 2.5
s_hat_4 = pinv(H)*y_4;
prior_mean4 = s_hat_4;
posterior_mean4 = prior_mean4 + lambda*(y_4_noise-H*prior_mean4);
figure()
imagesc(reshape(posterior_mean4, [8, 12])); colorbar;


