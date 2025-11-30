clc ; clear;
%%

% --- 1. Define Parameters ---
M_w = 4.3;   % Mass of one wheel
M_r = 10.12;   % Mass of robot body
m_b = 0.00271;  % Mass of the ball
R = 0.356;    % Radius of the wheel
h = 0.31665;     % Height of robot CoM from ball platform
r = 0.04006;   % Radius of the ball
g = 9.81;    % Gravity
d = 0.18865; 

% Safety limits
x_limit = d/4;                 % ball can move only within +/- d/2
theta_limit = deg2rad(15);     % 15 degrees in radians

% --- 2. Moments of Inertia ---
I_w = 0.5 * M_w * R^2;        % Inertia of one wheel (Note: Not used in M1, M2, M3)
I_r_com = (1/12) * M_r * (h^2 + d^2);       % body about its own COM
I_r = I_r_com + M_r * d^2;              % Inertia of robot body (parallel axis theorem)
I_b = 0.4 * m_b * r^2;      % Inertia of a solid sphere ball

% % --- 3. Derived Parameters ---
l = h + r; % Height of ball CoM from wheel center
% 
% 
M1 = [(I_b/(r^2)) + m_b, m_b*l;
      m_b*l,            I_r + m_b*(l^2)];


M2 = [0,             (g*m_b);
      (m_b*g), (M_r*d*g) + (m_b*g*l)];


M3 = -[m_b*l;
      I_r+ m_b*l^2];

% --- 5. Calculate A and B ---
% This is the numerically preferred method for A = inv(M1) * M2
A1 = M1 \ M2;

% This is the numerically preferred method for B = inv(M1) * M3
B1 = M1 \ M3;
%%
A = [0 1 0 0;... %[output:group:19492d60] %[output:8ea9a196]
    A1(1,1) 0 A1(1,2) 0;... %[output:8ea9a196]
    0 0 0 1;... %[output:8ea9a196]
    A1(2,1) 0 A1(2,2) 0] %[output:group:19492d60] %[output:8ea9a196]
B = -[0 ; B1(1); 0 ; B1(2)] %[output:88287b9c]
C = [0 0 1 0] %[output:758ffe6a]
D = [0] %[output:8c87c25a]
%%
% C_dash  = eye(4)
% D_dash = zeros(4,1)
% [num,den] = ss2tf(A,B,C_dash,D_dash)
[num,den] = ss2tf(A,B,C,D) %[output:63c095d8] %[output:09e92a9e]
%%
G1 = tf(num(1,:), den) %[output:7f3d3531]
%%
sys = ss(A, B, C, D);
%%
rank(ctrb(A,B)) %[output:6335ea38]
rank(obsv(A,C)) %[output:47b75a85]
%%
% step input
% step(sys)
%%
% sine input 
% t = 0:0.001:5;
% u = sin(2*pi*t);
% lsim(sys,u,t)
%%
% ramp input 
% t = 0:0.01:10;
% u = t;         
% lsim(sys, u, t)
%%
% LQR controller
% --- Bryson's-rule LQR weights (replace previous Q and R) ---
% State vector: x = [x; x_dot; theta; theta_dot]
x_max      = d/8;         % from your script
xdot_max   = 0.5;         % [m/s] <-- adjust if you expect faster/slower motion
theta_max  = deg2rad(15); % from your script
thetadot_max = 2;         % [rad/s] <-- adjust as needed
u_max = 5;                % [actuator units] pick actual motor torque/force limit

Q = diag([ 1/(x_max)^2,  1/(xdot_max)^2,  1/(theta_max)^2,  1/(thetadot_max)^2 ]);
R = 1/(u_max^2);

% Show numeric Q,R for confirmation
disp('Q = '); disp(Q); %[output:0b9c6ac1]
disp('R = '); disp(R); %[output:11e99c31]

% Compute LQR gain 
[K,S,P] = lqr(A, B, Q, R);

disp('LQR gain K ='); disp(K); %[output:36af7296]

%%
%K = place(A,B,[-8+5i;-9+20i;-9-20i;-8-5i])
%%
x0 = [0.01, 0, deg2rad(5),0]' %[output:1b06f103]
%%
% Closed-loop system
Acl = A - B*K;
sys_cl = ss(Acl, [], eye(4), []);

% Simulate response
t_sim = 0:0.001:20;
% Dynamics for ODE45
f = @(t, x) (A * x + B * (-K * x));
% Event function to stop simulation when limits exceeded
function [value, isterminal, direction] = stopEvents(t, x, x_limit, theta_limit)
    value = [...
        x_limit - abs(x(1));       % stops when this hits zero (ball position)
        theta_limit - abs(x(3))    % stops when this hits zero (tilt)
    ];
    isterminal = [1; 1];  % Stop integration for both events
    direction = [0; 0];   % Detect events regardless of sign
end

options = odeset("Events", @(t,x) stopEvents(t,x,x_limit,theta_limit));

[t, y, te, ye, ie] = ode45(f, t_sim, x0, options);
if ~isempty(te) %[output:group:6f0ae9ff]
    if ie(end) == 1
        fprintf('\n⚠️ Simulation stopped at t = %.3f s: |x| exceeded %.4f m.\n', te(end), x_limit);
    elseif ie(end) == 2
        fprintf('\n⚠️ Simulation stopped at t = %.3f s: |theta| exceeded 15 degrees.\n', te(end));
    end
else
    fprintf('Simulation completed without violating limits.\n'); %[output:84041628]
end %[output:group:6f0ae9ff]

figure;
plot(t, y);
xlabel('Time (s)');
ylabel('States [x, xdot, theta, thetadot]');
legend('x', 'xdot', 'theta', 'thetadot');
title('Closed-Loop Response with LQR (ODE45 + Safety Constraints)');
grid on;

%%
xb_max = max(y(:,1)) %[output:9791a366]
theta_max = max(y(:,3)) %[output:15857bb6]
% xb can be at max d/2 before the ball falls
disp_percent = (xb_max/(d/2))*100 %[output:6adab004]
rot_percent = (theta_max/(pi/2))*100 %[output:0e6468c6]
%%
%[text] ## Observer
G = 0.1*eye(4);      % or your actual process-noise input matrix
% Q = 0.2*eye(4);     % process noise covariance
Qdiag = 0.1*diag([10000000,100,1000000,10]);   % larger on states you want faster
Q_obs = Qdiag;
R_obs = 0.00001*eye(1);     % measurement noise covariance
N = zeros(4,1);               % no correlation
[L,P,E] = lqe(A, G, C, Q_obs, R_obs, N);   % or simply lqe(A,G,C,Q,R)
%steady‑state Kalman gain L, the error covariance P, and a vector E
%containing the eigenvalues of A-LC
Ke = L %[output:19a15f7d]
%%
x0_o = [0.025;0;deg2rad(5.25);0] %[output:973f237d]
%%
%[text] ## Min order observer
% Measured and unmeasured state indices
ia = [3];  % measured: theta, thetadot
ib = [1, 2, 4];  % unmeasured: x, xdot

% A partitions
A_aa = A(ia, ia);
A_ab = A(ia, ib);
A_ba = A(ib, ia);
A_bb = A(ib, ib);

% B partitions
B_a = B(ia);
B_b = B(ib);

% C partitions
C_a = C(:, ia);
C_b = C(:, ib);

% =========== DYNAMIC DIMENSIONING USING ia, ib ===========
na = length(ia);  % number of measured states
nb = length(ib);  % number of unmeasured states

%%

% =========== DYNAMICALLY SIZED MATRICES ===========
G_min = 0.01 * eye(nb);    % nb × nb  process noise input
Q_min = 10000 * eye(nb);     % nb × nb  process noise covariance
Q_min(2,2) = 100 %[output:63b81036]
R_min = 0.005 * eye(na);   % na × na  measurement noise covariance
N_min = zeros(nb, na);     % nb × na  cross-covariance (typically zero)

% Design reduced-order Kalman filter
[L_min, P_min, E_min] = lqe(A_bb, G_min, A_ab, Q_min, R_min, N_min);

Ke_min = L_min %[output:9d11ed4b]
%%

A_hat = A_bb - Ke_min*A_ab;
B_hat = A_hat*Ke_min + A_ba - Ke_min*A_aa;
F_hat = B_b - Ke_min*B_a;
na = length(ia);   % measured states
nb = length(ib);   % unmeasured states

C_hat = zeros(4, nb);   % from observer state z to x_hat
C_hat(ib, :) = eye(nb);

D_hat = zeros(4, na);   % from measurement y to x_hat
D_hat(ia, :) = eye(na);

% Check observability
obs_rank = rank(obsv(A_bb, A_ab)) %[output:0ea290d6]
%%
x0_omin = [-0.02;0;0] %[output:3c41c227]
%%
%[text] ## FFT
% y_noise  = out.y_noise.signals.values;
% ti = out.y_noise.time;
% save('noise_data.mat','y_noise','ti');
%%
data = load('noise_data.mat') %[output:1f31598f]
y_noise = data.y_noise;
ti      = data.ti;
Ts = ti(2) - ti(1) %[output:6c5a19e0]
Fs = 1/Ts; 
N = length(y_noise);
Y = fft(y_noise - mean(y_noise));        
f = (0:N-1)*(Fs/N);          
P = abs(Y).^2 / N;           
figure; 
plot(f(1:N/2), P(1:N/2)); %[output:91924c23] %[output:49a842be]
xlabel('Frequency (Hz)');
ylabel('Power');
grid on;
title('Noise FFT');

%%
%[text] ## Low pass filter
fc = 20;
tau = 1/(2*pi*fc) %[output:22ddab3e]
num = [1] %[output:0f6679e5]
den = [tau   1] %[output:97a12f32]

% 10 Hz is far above the 1.5 Hz interference
% the robot's physical angle dynamics are active around 3–6 Hz → fully preserved
% it significantly reduces high-frequency measurement noise
% it barely introduces delay (a few milliseconds)
%%
%[text] ## Butterworth filter
fc  = 20;        % cutoff
Fs  = 1/Ts;     % sample rate from your model
Wn  = fc/(Fs/2);
n   = 2;        % 2nd order

[b,a] = butter(n, Wn, 'low');


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:8ea9a196]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["-0.0143","0","-3.0441","0"],["0","0","0","1.0000"],["0.0560","0","39.4486","0"]]}}
%---
%[output:88287b9c]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["0"],["0"],["1"]]}}
%---
%[output:758ffe6a]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":1,"type":"double","value":[["0","0","1","0"]]}}
%---
%[output:8c87c25a]
%   data: {"dataType":"textualVariable","outputData":{"name":"D","value":"0"}}
%---
%[output:63c095d8]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":1,"type":"double","value":[["0","0","1.0000","0","0.0143"]]}}
%---
%[output:09e92a9e]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"den","rows":1,"type":"double","value":[["1.0000","0.0000","-39.4343","-0.0000","-0.3923"]]}}
%---
%[output:7f3d3531]
%   data: {"dataType":"text","outputData":{"text":"\nG1 =\n \n                      s^2 + 0.01427\n  ------------------------------------------------------\n  s^4 + 1.905e-15 s^3 - 39.43 s^2 - 5.046e-15 s - 0.3923\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 49 32 48 32 48 46 48 49 52 51 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 49 46 57 48 52 55 101 45 49 53 32 45 51 57 46 52 51 52 51 32 45 53 46 48 52 53 56 101 45 49 53 32 45 48 46 51 57 50 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:6335ea38]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:47b75a85]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:0b9c6ac1]
%   data: {"dataType":"text","outputData":{"text":"Q = \n   1.0e+03 *\n\n    1.7983         0         0         0\n         0    0.0040         0         0\n         0         0    0.0146         0\n         0         0         0    0.0003\n\n","truncated":false}}
%---
%[output:11e99c31]
%   data: {"dataType":"text","outputData":{"text":"R = \n    0.0400\n\n","truncated":false}}
%---
%[output:36af7296]
%   data: {"dataType":"text","outputData":{"text":"LQR gain K =\n -211.3992 -131.3312  162.7202   18.2124\n\n","truncated":false}}
%---
%[output:1b06f103]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0","rows":4,"type":"double","value":[["0.0100"],["0"],["0.0873"],["0"]]}}
%---
%[output:84041628]
%   data: {"dataType":"text","outputData":{"text":"Simulation completed without violating limits.\n","truncated":false}}
%---
%[output:9791a366]
%   data: {"dataType":"textualVariable","outputData":{"name":"xb_max","value":"0.0100"}}
%---
%[output:15857bb6]
%   data: {"dataType":"textualVariable","outputData":{"name":"theta_max","value":"0.0873"}}
%---
%[output:6adab004]
%   data: {"dataType":"textualVariable","outputData":{"name":"disp_percent","value":"10.6016"}}
%---
%[output:0e6468c6]
%   data: {"dataType":"textualVariable","outputData":{"name":"rot_percent","value":"5.5556"}}
%---
%[output:19a15f7d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"4","name":"Ke","rows":4,"type":"double","value":[["2.9501"],["-0.1367"],["1.0001"],["0.5787"]]}}
%---
%[output:973f237d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_o","rows":4,"type":"double","value":[["0.0250"],["0"],["0.0916"],["0"]]}}
%---
%[output:63b81036]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"Q_min","rows":3,"type":"double","value":[["10000","0","0"],["0","100","0"],["0","0","10000"]]}}
%---
%[output:9d11ed4b]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"Ke_min","rows":3,"type":"double","value":[["18.3261"],["0.2491"],["14.2145"]]}}
%---
%[output:0ea290d6]
%   data: {"dataType":"textualVariable","outputData":{"name":"obs_rank","value":"3"}}
%---
%[output:3c41c227]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_omin","rows":3,"type":"double","value":[["-0.0200"],["0"],["0"]]}}
%---
%[output:1f31598f]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"data","value":"         ti: [100001×1 double]\n    y_noise: [100001×1 double]\n"}}
%---
%[output:6c5a19e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"Ts","value":"1.0000e-04"}}
%---
%[output:91924c23]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Integer operands are required for colon operator when used as index."}}
%---
%[output:49a842be]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Integer operands are required for colon operator when used as index."}}
%---
%[output:22ddab3e]
%   data: {"dataType":"textualVariable","outputData":{"name":"tau","value":"0.0080"}}
%---
%[output:0f6679e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"num","value":"1"}}
%---
%[output:97a12f32]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"den","rows":1,"type":"double","value":[["0.0080","1.0000"]]}}
%---
