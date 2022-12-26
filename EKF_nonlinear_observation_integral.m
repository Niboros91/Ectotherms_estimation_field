function [xhatOut, P] = EKF_nonlinear_observation_integral(iteration,A,L,C,x_ini,y,measuring,P)
% Define storage for the variables that need to persist
% between time periods.
%measruing  -> 0,1 1 if measurement
% y -> measurement value

persistent xhat Q R 
 

if iteration==1 %In the first iteration we create the 

    xhat=x_ini;
    R = 0.1;     % Only one state
    Q = eye(16)*0.02;
   
end

% Estimation Kalman filter
xhat = A*xhat;
P = A*P*transpose(A) + L*Q*transpose(L);
P=(P+P.')/2; %To force symmetry
P_pre=P;


if measuring == 1
    % Update Kalman Filter
    H = xhat(9);
    K =(P*C')/(C*P*C' + H*R*H'); 
    % K =(P*C')/(C*P*C' + R); 
    xhat=xhat+K*(y-C*xhat);
    P=P-K*C*P;
end

% Post the results
xhatOut = xhat;


if measuring == 1
    % Update Kalman Filter
    xhat(9) = 0; %If measurement we move id to 0
end
end