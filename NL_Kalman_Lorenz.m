clear, close all
%% PART I:  Understanding Lorenz

randn('seed',123)

% For these given parameters the system exhbits chaotic behavior
sig = 10; rho = 28; b = 8/3;

% lorenz system is a system of 3 ordinary differential equations
% Y(1) = x, Y(2) = y, Y(3) = z
Lorenz_equation = @(t,Y,sig,rho,b) ...
    [sig*(Y(2)-Y(1)); rho*Y(1)-Y(2)-Y(1)*Y(3); Y(1)*Y(2)-b*Y(3)];

Lorenz_Jacobian = @(Y, sig, rho, b) ... 
    [-sig, sig, 0; rho-Y(3), -1, -Y(1); Y(2), Y(1), -b];

D = b*rho-b;
% if D < 0 there is only one equilibrium point in the origin
if D<0
    ne = 1;
    xe = 0; ye = xe; ze = xe^2/b;
    disp('There is one simple equilibrium point:')
    disp(['x = ',num2str(xe),'  y = ',num2str(ye), '  z = ', num2str(ze)])
% a pitchfork bifurcation occurs at D =0
elseif D==0
    ne = 1;
    xe = 0; ye = xe; ze = xe^2/b;
    disp('There is one triple equilibrium point:')
    disp(['x = ',num2str(xe),'  y = ',num2str(ye), '  z = ', num2str(ze)])
% 2 additional critical points appear
else D>0
    ne = 3;
    xe = [0;sqrt(D);-sqrt(D)]; ye = xe; ze = xe.^2/b;
    disp('There are three equilibrium points:')
    disp(['x = ',num2str(xe(1)),'  y = ',num2str(ye(1)), '  z = ', num2str(ze(1))])
    disp(['x = ',num2str(xe(2)),'  y = ',num2str(ye(2)), '  z = ', num2str(ze(2))])
    disp(['x = ',num2str(xe(3)),'  y = ',num2str(ye(3)), '  z = ', num2str(ze(3))])
end

% The eigenvalues of the Jacobian at each equilibrium point determine
% whether the equilibrium is stable or not. The real parts of all
% eigenvalues have to be negatives in order to be stable.
for k = 1:ne
    J = Lorenz_Jacobian([xe(k), ye(k), ze(k)], sig, rho, b);
    Le = eig(J);
    if (real(Le(1))<0)&(real(Le(2))<0)&(real(Le(3))<0)
        disp(['Equilibrium at point ',num2str(k),' is stable.'])
    elseif real(Le(1))*real(Le(2))*real(Le(3))==0
        disp(['Equilibrium at point ',num2str(k),' is (neutrally) unstable.'])
    else
        disp(['Equilibrium at point ',num2str(k),' is unstable.'])
    end
end

% Take an initial guess
x0  = -3.3; y0 = 3.0; z0 = 10.5;

time = 4.0;

options = odeset('RelTol', 1.E-7, 'AbsTol', 1.E-7);
[T,Y] = ode45(Lorenz_equation, [0:0.002:time], [x0;y0;z0], options, sig,rho,b);


% plot the state space of lorenz equation
figure(1), clf
comet3(Y(:,1), Y(:,2), Y(:,3)) % trace curve through the points [xi,yi,zi]
plot3(Y(:,1), Y(:,2), Y(:,3),'DisplayName','state-space'), grid on
xlabel('x'), ylabel('y'), zlabel('z');title('Trace Curve of Equilibrium Points');
legend('show');
disp(['at time  ', num2str(T(length(T))), '   the final point is:  '])
disp(Y(length(Y),:))


figure(2), clf
subplot(3,1,1)
%plot(T,Y(:,1))
subplot(3,1,2)
%plot(T,Y(:,2))
subplot(3,1,3)
%plot(T,Y(:,3))

%% PART II:  Lorenz with noise
%generate data for "actual" system with noise
Dt = 0.01;
Q = 0.01*eye(3); %small variance
% Q = 1*eye(3); %large variance
B = chol(Q)'; 
k_steps = 200;
Xt = zeros(3,k_steps);
Xt0 = [x0; y0; z0]; % we will use the same initial point as the deterministic simulation

% in this part we will simulate the lorenz equations while adding
% white noise from time 0 to time 2 at a time interval of 0.01
for k = 1:k_steps
    Sol = ode45(Lorenz_equation, [0 Dt], Xt0, ...
        options, sig,rho,b);
    Xt(:,k) = deval(Sol, Dt) + B*randn(3,1);
    Xt0 = Xt(:,k);
end
% Plot the "actual" data (state + deterministic + noise)
figure(1), hold on
comet3(Xt(1,:), Xt(2,:), Xt(3,:))
hold off

figure(2)
subplot(3,1,1)
hold on
plot([Dt:Dt:k_steps*Dt],Xt(1,:),'r','DisplayName','actual equilibrium x');
hold off
subplot(3,1,2)
hold on
plot([Dt:Dt:k_steps*Dt],Xt(2,:),'r','DisplayName','actual equilibrium y');
hold off
subplot(3,1,3)
hold on
plot([Dt:Dt:k_steps*Dt],Xt(3,:),'r','DisplayName','actual equilibrium z');
hold off

% generate data/observations, only the x coordinates are observed
% sqrt(r)*randn term provide the noise to the system
H = [1 0 0]; R = 1;
for k = 1:k_steps
    y_data(k) = H*Xt(:,k) + sqrt(R)*randn; 
end
    

%% Part III. Extended Kalman filter.
% Invert the observations in order to solve for equilibrium points

Q = 10000*Q; 
X_hat = zeros(3,k_steps);
res = zeros(1,k_steps);
sig2 = zeros(1,k_steps);
stableq = zeros(1,k_steps);
X_hat0 = [x0; y0; z0];  % we will still use the same initial point
SIG0 = zeros(3); 

for k = 1:k_steps
    disp(num2str(k))
    %prediction
    Sol = ode45(Lorenz_equation, [0 Dt], X_hat0, ...
        options, sig,rho,b);
    X_hat(:,k) = deval(Sol, Dt);
    
    A = Lorenz_Jacobian(X_hat0, sig, rho, b);
        %stability of underlying dynamic system, if the real part of
        %eigenvalue of lorenz Jacobian is negative, the nthe system is stable
        d = eig(A);
        if max(real(d))<0
            stableq(k) = 1;
        end
   
    PHI = expm(A*Dt);   %PHI is the approximated transition matrix
    SIG = PHI*SIG0*PHI' + Q; % specifically covariance sig(k+1|k)
    
    % update
    sig2(:,k) = H*SIG*H'+R;     % residue variance
    K = SIG*H'*inv(sig2(:,k));  % Kalman gain
    SIG = SIG - K*H*SIG;        % specifically SIG(k+1|k+1), posterior covariance
    res(:,k) = y_data(k)-H*X_hat(:,k);  % residue = y - H*x_hat
    X_hat(:,k) = X_hat(:,k) + K*res(:,k); % prediction after Kalman correction
    
    X_hat0 = X_hat(:,k);
    SIG0 = SIG; 
        %%stability of feedback system           
%             d = eig(PHI-K*H*PHI);
%             if max(real(d))<1
%                 stableq(k) = 1;
%             end
end

% In figure 2 we plot the prediction vs. actual
figure(2)
subplot(3,1,1)
hold on
plot([Dt:Dt:k_steps*Dt],X_hat(1,:),'k','DisplayName','predicted equilibrium x')
xlabel('time'); ylabel('x coordinate');hold off; legend('show');
title('actual equilibrium x vs. predicted equilibrium x')

subplot(3,1,2)
hold on
plot([Dt:Dt:k_steps*Dt],X_hat(2,:),'k','DisplayName','predicted equilibrium y')
xlabel('time'); ylabel('y coordinate');hold off; legend('show');
title('actual equilibrium y vs. predicted equilibrium y')

subplot(3,1,3)
hold on
plot([Dt:Dt:k_steps*Dt],X_hat(3,:),'k','DisplayName','predicted equilibrium z')
xlabel('time'); ylabel('z coordinate'); hold off; legend('show');
title('actual equilibrium y vs. predicted equilibrium y')

% figure 3 plots the normalized residue and the underlying system stability
figure(3), clf
subplot(2,1,1)
plot(1:k_steps, res./sqrt(sig2))
xlabel('time'); ylabel('normalized residual'); title('normalized residue');
subplot(2,1,2)
area(stableq)
colormap summer
title('green means stable')
xlabel('time'); ylabel('stability (1 is stable)');

% figure 4 plots the residual variance and underlying system stability
figure(4), clf
subplot(2,1,1)
plot(1:k_steps, sig2)
xlabel('time');ylabel('residual variance'); title('residue variance');
subplot(2,1,2)
area(stableq)
colormap summer
title('green means stable');
xlabel('time');

% This part performs the post-factum check on the estimator.
figure(6), clf
plot(1:k_steps, res,'r.')
title('residuals'), xlabel('time'), ylabel('residual');

mres = mean(res)
sdres = std(res)
figure(7), hist(res);
xlabel('reisdue magnitude');;title('residue histogram')
corrcoef([res(2:end)', res(1:end-1)'])
