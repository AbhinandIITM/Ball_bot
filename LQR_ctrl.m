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

figure; %[output:2b3d5351]
plot(t, y); %[output:2b3d5351]
xlabel('Time (s)'); %[output:2b3d5351]
ylabel('States [x, xdot, theta, thetadot]'); %[output:2b3d5351]
legend('x', 'xdot', 'theta', 'thetadot'); %[output:2b3d5351]
title('Closed-Loop Response with LQR (ODE45 + Safety Constraints)'); %[output:2b3d5351]
grid on; %[output:2b3d5351]
%%
function info = impulseinfo(t, y)
    % simplified impulse-response info
    % fields: Peak, SettlingTime, Overshoot

    y = y(:);
    t = t(:);

    % Peak
    [Peak, ~] = max(abs(y));

    % Final value
    final_val = y(end);

    % Overshoot (only meaningful if final != 0)
    if abs(final_val) > 1e-6
        Overshoot = (Peak - abs(final_val)) / abs(final_val);
    else
        Overshoot = NaN;   % consistent with control theory practice
    end

    % Settling time (2% band)
    err = abs(y - final_val);
    tol = 0.02 * abs(final_val);
    idx_set = find(err <= tol, 1, 'first');
    if ~isempty(idx_set)
        SettlingTime = t(idx_set);
    else
        SettlingTime = NaN;
    end

    info = struct('Peak', Peak, ...
                  'SettlingTime', SettlingTime, ...
                  'Overshoot', Overshoot);
end
%%
impulseinfo(t,y) %[output:9791a366]
%%
xb_max = max(y(:,1)) %[output:15857bb6]
theta_max = max(y(:,3)) %[output:6adab004]
% xb can be at max d/2 before the ball falls
disp_percent = (xb_max/(d/2))*100 %[output:0e6468c6]
rot_percent = (theta_max/(pi/2))*100 %[output:195c8ce3]
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
figure;  %[output:8be810fa]
plot(f(1:N/2), P(1:N/2)); %[output:91924c23] %[output:49a842be] %[output:8be810fa]
xlabel('Frequency (Hz)'); %[output:8be810fa]
ylabel('Power'); %[output:8be810fa]
grid on; %[output:8be810fa]
title('Noise FFT'); %[output:8be810fa]

%%
%[text] ## Low pass filter
fc = 14;
tau = 1/(2*pi*fc) %[output:22ddab3e]
num = [1] %[output:0f6679e5]
den = [tau   1] %[output:97a12f32]

% 10 Hz is far above the 1.5 Hz interference
% the robot's physical angle dynamics are active around 3–6 Hz → fully preserved
% it significantly reduces high-frequency measurement noise
% it barely introduces delay (a few milliseconds)
%%
%[text] ## Butterworth filter
fc  = 14;        % cutoff
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
%[output:2b3d5351]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAe4AAAEqCAYAAADTbrP0AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ20VsV57x8VkeNHlIMKwuFDLTGtLZg0LeQkXFMT+5ErpMYkB44Wi4SyevFgGr5PY4w3KSDiapCYW6LkqEsFmjRUMG1T2zSB5FR6lYCJqYYb5eOInCCiUYOgkbuerfNmzpzZe8\/M3vud\/e7577VcAu\/MPDO\/Z\/b893yfcPz48eOEBwRAAARAAARAoCEInADhbgg\/IZMgAAIgAAIgEBGAcKMigAAIgAAIgEADEYBwN5CzkFUQAAEQAAEQgHCjDoAACIAACIBAAxE44fzzz48Wpy1btoza2tr6ZP2WW26hNWvW0OzZs2nRokWlK9YLL7xAM2fOpJ07d9IDDzxAEydOjM3jrl27aMaMGdHvXV1dNHbs2NKVJy5Dcjnjwuj81zAFbNCMbtiwgZYsWUKTJ0+m5cuXU1NTU1QS9te8efOos7Mzqmc29VSHwqXuxtUZ3XuSZ\/0SeWUmapvxyCOPUHt7e58iDh8+XPs+CrZxVUMth0kZOK24dkLkbfz48bR27Vpqbm6umY7LSz3aRdW2Ln+mr8+RI0do8eLFtHnz5ihKHPu49ET8adOmJba1pvkR4USd2b9\/fy1qlnLa2jcNX0T5Xd5tzm9NuPkvKiwIt6lLiw1n2ihBvIv1g5q6Trh1L2K9hVu8t3E0dB8a4gM4iaBJ\/RK2VYFMy5MqgmnCzfmU45i+I2kfLjrBiMt70cIdZ9dWcIVPs3wEZK3DcfUqrV6Y1Ll6vPVFld9ZuHlVuexQuTJCuOtRJdJtyJUmaWSkjF+p6aWrVgjfwi33apPEUK5HedUvUfZzzjmnT69Vbl\/U+is33PJvcaMZXFvkOEKIszSscnrqOyR6WY899pjzSJ3wiY3g6nwi\/5vLR0MS07S3MAvfuLTluhpXH22YpZUhy+9FlD9Lfmpz3Fx5+ZGHt+KEW\/d1qxMU3Vdz2hevruev+2Jkh65cuTJ6iYsYKjctI+dNN9QT1+PgF27MmDHRECs\/JhUzrWGNG+ZTy6ATdl051R6ZLEYyc7XXI1dE9Us6qZfHrL73ve9F0zJxTHSMdY2XSZnlfMa9kCL\/MjNVlLdv395nqJzroToUzOVeuHAhXX\/99bV6umfPHiv\/m36Vy0OhcQ27eCflcrnWr7gRiLiecFzvScc6SWR05XRtWNXhe\/UdEelyWdUhdNOG10W45foutyVJ0zPqqInwgzpELvIt+0PlIPswbiib6\/YXv\/jFfiO1sr04n6fVVbbJ782KFSv6TWumjRyodSFr25I0lH\/nnXfWppOZK7dhcpuu5lVt73Xvtuzj6667jubMmUNiGkH4JXFxmk64dfNUoiLIzo4b6orLuDy\/oWu804bO8pzjNi0j5zMpX3KlTRsSSsp\/WsMq0pbFUVfZVK5Jw4tpL67caMl2k9KMEwtdA6gTTLWOcDyZsUmZdbbUei6XQa6vol6I8m7atMlauCdNmkRbt27tk420kRJT4Y5r7NM+VFzql8oxrrE2Ea20xkteP6B+xAt2\/O+m611EGnJvmuN+4Qtf6CdCcW1BWnsj8zFhkMQzqTPDv6W9x3Pnzu0zt60Kd1wbJup5T09PtD5InYO+9dZbacGCBf06Tib11aSupjGRf5ff07RpE9u2xUS4de2haB\/Ucsh5Tar7cR+G3OZZCbf8csqNtVy5RYVWG0PdF1bcyy4qkhCPuCEi3ZBZXGFNKhPHtSljXL5E\/mUHpQ3vqT3SuMY26StfsDflqmOia2TkiivnU37hhW1d2XWc5H+TXyRdXdKNKKg9R14YJhbe6IZck4YWVUFWGyqRnvCh+LvLHHdaWXX+Na27JuF0dSOtoRN5Mv24lMMlLfoS6ep6y2nDumq6snDHvSPqOya3M5deemk0WqJ+RCV9mCfxSOtocB7TxD\/uQ1SNp7aXnLbuY960vsp1JG0qQte50+VH9YlJXVXj6NoGXXvN\/yY+4tLeN5O2hRcqpo3Mqf4U4Q8ePFibYtF9rKQJt2hr1HJaCXccbJ0oqxXXZFhTdZSA\/tOf\/jR6qdTeus3wmGlFsSmjyRC12uirHGwbNl2jpKaZ1hDreimcbtywvQkTLueUKVO04slpq0IrN7Rx81ui0VAbMNMhf5lVUq9WfcnEELjoHTPfWbNmRY2B\/CKaNoRxPVrT+pu17socXITbZD41Lo+29Vv9AIz7oM0q3CK\/Yj5etDFxC3R1o2dJdSoP4RZ+043WpU2FiLi6D+24TldS2xJXV5P8kPRhYlqn5TzFTd2qnY0hQ4bUhDuPtkUd1ZDLpftAUjnq2mO1beM4YreTbjpLbkPZf1bCnfQS6qDGDQ8LmHFflKLgqnDHzT\/Jc9yqTRHn0KFDRtvBbMoY1yvQNY5plc5EWLicOnZqXFOu\/CUZFzZpbldspVPLefnll8cOV5q+4GkNhPpSJA3nmQq33EO58cYbaceOHdGWmS9\/+cvEc1j88NqP+fPnk7zwykW45Ze+KOHmIc24RlP3EaH7N7lepA3lM5804U5ay1HvoXJdj9LkA0PUJ9uhb9vwcSMGcr51w626qaQ04U77wNCNesp1S2XJeedOVtIIolxfkupq2gen+F2tP7Jwm7xvcdMhcVOApsKdND1qItxxI0R9hJuNMER5PkkVG5OeV9yXoAzHdl4qzq6u4csq3DZldOlxq84waTDiemwyU5NKFtcgiH9XF7HoPhLiXtosPW6Tl0vOu\/oFy\/Hf+c53Ws9xymnKfuBeNQs0CzcvjmERl3vfYgFn2YQ7bsEPvxNcPl5c9fDDD9cWxumGQeMWLJk2xHLPQe2p5LU4Tf7QShMWXZ1P+7DlOJxX+UNUx8VkYSmn5SLccb0utc3jtlRMEckfWKZD5SbtT1KPU+4Jsi\/44QVaadu4TBan8bw6L4YV53OkdX6EP2yFO61tYftpQ+VxI6miLvEZKbZD5anCLX91JW0HkzOfNMctVyY5PbUCt7S01CqdnF7SYiFd\/hiO6XyR2rCoL7ZpGdmZsiN0C\/OKmONWXwjZd3HzIeJjTOWqa1B0L5TNHLewEbdYJK2hTZrvTFoJLQ\/TJ9WluI8X9WNA5FPtkaSt8NV9+MW99Hn3uGWR4D\/rWIvyx33omdQvHcOkssS1L7IAyw2cLAa6Dwbd2gpTlmpvL64+qFNcaaKY9FHsItzyO6dbs6ETKREurg1L+tDkTpsuftoct46nyQiNWlfjhrTjVmjr1vOIuiLPcad1CnQfSGmjUroed9zWS5mFrt7ajDbJ\/rM+gCVuWEFuKFSnqJVaJ+ZqGPVrNm1Ix1S4dcNJwrZIw7SMcgOje3FNV5WnfbnH9bjZpiy0cQsx5LzJtuK2iXB406E41e9J8+txwpv2cpmmGee3NL6CjzxaYzIlkdYQcrpc5rjVt6ZiY9JD1DUmSWISx1wV7rj6paadtgUobVdFXMOXVIa4bWdJceJ6gnE9zzj2pnUqKS9pvyUxi+vAqGmmDZWrH09y\/KTdImr55byarInQvXNpbahpe2XzoWzatqjhRPm\/+c1vak8XTXtnMw+ViyNPdbDjhiaSJttl+LrM6+yo4UyOQXTdx20i3FwG0zLqvjh1+ZdZihWsqkDGvchJwq3aT9rGFcdVbSCS9hrecccd9LWvfa12bGLcS6qmmbSALk241Q8UwUn3ZW9al3SsdSvik4b04tY4yB+ZPoQ7rv6aNMo6YYsbWVIZpq0m1rUHJu+6akcXJ6kBluPbCreOZT1EW+RZ12FJ+rgSH9If+9jH+q3piauvHEe1oxvpkD+MVQa6Fd9pHybid5t6YZJXG+G2aVt05Y8Tbk7XJLxzj7us93Grvae0HrVcSURjm\/eZuqYVURcu7iMoS5r1iuuyArReeYMdewLCnx0dHf3uJ7BP7dcx1JXa8nnfWdJF3PITSPooKH\/uGy+HpbwdjEWbV\/CK5fHq39MwC5G0Efu0NLP+DuHOShDxG4FAGd+9RuDWqHlUR9bSFqU1ajnLlu\/SCbfoLfNQjHr8KsNLuqVMt9o46cawejoDwl1P2rDli4DodfNQaxlvFPTFpap2Xee2q8qjXuUqnXAL8eWXXhZd7nVzJYk7L1jEGzVqFInzXeWtBPUCCjsgAAIgAAIgUCSB0gk3f7HrDpe3GS4XX\/0Q7iKrDtIGARAAARDwQcBIuE1XbOoKYLqnT8SFcPuoBrAJAiAAAiDQKASshFsdvk4rZNrwti5+vYWbL5TAAwIgAAIgEA4BPvyrkZ+ghZtFm6+l27ZtWyP7EHkHARAAARCwIDBhwoToYKRGFXAj4bbgkTmo6+I02bDpHLfYK84OHDFiROa8IwFzAvyxtGrVqujlAXtzbnmFBP+8SLqlA\/5u3PKIJdhv2bIFwp0HUE4jy3YwkQdb4W5kB+bFvd7piI8msK83+bfsgb8f7sIq+PvjXwX21j1u7hHPmzePOjs7SVzvKLuAT9BZv3597LYtE3cJsPLZ4fKBLGlpQLjTCPn\/fffu3XT33XfTZz\/7WRowYID\/DAWWA\/D363Dw98cfwj12bD\/6eQi33CMQBtRT0HhPNz+6Qx5shHve9e204R8bd8jEX\/XPZvm1116j5557jkaOHAnhzobSKTb4O2HLLRL454bSOqFghDvpVhYdtUY59u715x+hl7un0Vkf3konntrYqwyta6\/nCGi4\/DoA\/MHfLwF\/1oMRbhlx2lC5P3fYW2YHvvPn0+i0d99Kp4z8uH0CiOFMAMLhjC6XiOCfC0bnRMDfGV3miEEKd2ZqJUqAHXj6k9Po\/PEfj8Tb9Hlp37Eo6JkjB5pGQTiFABouv1UC\/MHfLwF\/1oMWbts7VP25Kd4yO\/Cph6bRJ\/\/sBmq66NNGWWTRfnz9YdrT\/QqNbj2dJi0YahQPgfoSgHD4rRHgD\/5+CfizHqxwi4Krc9niTtYyXaeZVD2EcH\/0gy3RPLfJw4K99dbeSLQf33CYrri9JfozHjsCEA47XnmHDpV\/WU5KfOONN6LFmeeddx4WZ2as3LaHqAQp3HH7rAV7Xu29f\/9+Wr58OTU1NWV0SbHR2YEPfKmdPt9+nJqnPGNkjEWbxfuajRdGAr51ZS919o4ziotAvyYQqnCUpQ6EyB8nJZal9uWbD9tT0IIUbnGy2dSpU6mtra2fB\/LaDpava\/WpycJtsrKch8nvu\/LpPr3spUMfpytuH0nj2gbXI8uVsRGicJTJeSHyx0mJZaqB+eTF5RS0IIW7aj1u3se9+XPH6YzWdXTy2RMTaxP3tO+\/8mm6euMFteHx+678GZ01cmAk3njMCYQoHOZ0ig8ZIv8qNNjF14zGsuDiU5c4ZaNifXIaF6BKc9zt7e306JeOG20J46Fxntee8+i7an5kMX9obg9ds\/ECrDK3qN0hCocFnsKDhsi\/Cg124RWjwQy4+NQlTtmwOAk3F6Iqq8pZuDd97jiN\/VD6Xm7uXfPD89vyc8d7n4xWl2O43Lx6hygc5nSKDxki\/yo02MXXjMay4OJTlzhlo+Is3GUriEt+hAO\/tXwEjRo7MXUvNwu3bgvYQ3P3Rb1tbA0z90KIwmFOp\/iQIfKvQoNdfM1oLAsuPnWJUzYqEO72dvr2nVfR0EE99I73r4v1Dy9M4561biEaD5\/zMDqGy82rd4jCYU6n+JAh8q9Cg118zWgsCy4+dYlTNioQ7reF+5zjjyTu5RbCLS9Mk52J1eV2VTtE4bAjVGzoEPlXocEutlY0XuouPnWJUzYyRsIttoDt3LnTKP\/jx4\/PdK2nkZEcAgkHPnz\/DTT44JcS93KLFeW8ME131ClWl9s5JEThsCNUbOgQ+VehwbapFbw1d\/Xq1dTV1RVdwZy2I8gm7bKEdfGpS5yylFfkw0i41UyLhWkdHR199nKrFaVshVXzowp30l5uMRwuryiX0+Pfea4bh7GYeT1E4TAjU59QIfKPa7D3vvBafaAXbGVU86B+FuTrj7l97u7ubojDsUxRuYiwSxzT\/NQrnLVwp321NdrJabyq\/Lv\/8kB02UiScOu2gqlO4uHyuB55vRzaKHZCFI4y+SZE\/nEN9i3ffoZu+fbuMrnHKS+L\/mgMLfqj8\/vEFZ2smTNn0qZNm2jFihVR77sqj4sIu8QpGy9r4a7ayWks3Fu2bKFTt09K3MsdtxVMdiiGy82rd4jCYU6n+JAh8k\/qcVeh1809bl2vW9whod4tUXwtK96Ciwi7xCm+JHYWrIW7ij1uFm7ucTdddEPsvdwmosxD5TxkjuHy9EoYonCkU6lfiBD5V6HBdqkhQrhnz55NixYtckmitHFcfOoSp2wArIWbCyAKrrsdTF4MUbbCqvmRHcjCfcrIq2Kv9+StYHzAStJebTHPjeHydM+HKBzpVOoXIkT+VWiwbWsID5UvXLiQZs2aRUuXLqWVK1fSxInJRzvb2vAZ3sWnLnF8llFn20m4OSHdSvPJkyc31MIH2YGDD66K+Jz27lu1PjLd7vXWXm9c9ZlW0UMUjjQm9fw9RP5VaLBt6ogYHW1tbY0WETfSBVCm5XTxqUsc0\/zUK5yzcNcrg0XaUYX7V7\/UH8KStodbziMvYuPwuHQk2XMhCkeRddk27RD5V6HBtvEzLxTmMq9du5aam5uxHexteFWoBxDutxenDXn1G3TkqVXavdw2wp2239vmxaty2BCFo0z+DJF\/FRrsMtWhMuTFxacuccpQVjkPEO63hZtPTnv1hwu0wm0rxjysHnfCWtkqgK\/8hCgcvljr7IbIvwoNdpnqUBny4uJTlzhlKGtm4U47Sa3RTk7jVeV8VvnL3fq93LaHq\/Dqcn4wXB5f3UMUjjK9\/CHyr0KDXaY6VIa8uPjUJU4ZyppZuOVDVnhT\/+7du6NtBmIFY6Ns8lcd+MKm87V7uU0OX5GhYnV5ejUPUTjSqdQvRIj8q9Bg16+GNIYlF5+6xCkbDeuhcvUAFoawbt262mpyXrkohLxshVXzU5Rwsx3TVehlZ1RU\/kIUjqJYuqQbIv8qNNguvq5yHBefusQpG8PMws29bN4feNttt0UrF9W\/17PA6hB+2oEDqgNf\/LdJdPKQ\/vdy89D3i\/uO0TUbLzQuDsfZ0\/1qdAQqnv4EQhSOMtWDEPlXocEuUx0qQ15cfOoSpwxlzTRUru4NZLGcN28edXZ2Rmfg+hJudSQg7WhWhqA6kBen8aPu5TY57lTnWOzpjq\/uIQpHmV7+EPlXocEusg7ZTHUKHZg2bZrXA11cfOoSp0juLmlb97jZiLqRn+e8x4wZU9vk7+MGGt3hAuwgzpvYx5g2VM7CrdvLbXLcKaf9xsG3LioYcM6Y6P8cb3Tr6Ymnrbk4rQpxQhSOMvktRP5VaLCLrEMQ7iLp5pu2k3BzFuQFavz1xbfP8H3dvlaUy9fXCUSi180L53TH\/PGLPOUrO2jChAnRrTq\/d8q\/abeEpR13euSJ79Lhr98cmX3j52+Jd9PFH6RneqfQT75zJl257gPaO7zzdWVjpRaicJTJQyHyD024Te7j5nZzzZo1UdXkC5eeeOKJ2g1iole9efPm6Hcx9aj+u8\/LS1x86hKnTO8u58VZuMtUEHX4XhXuqVOn9rk3XPzODvzk3JupY8kXaNnWl+nTv\/UYXX\/uHf2u94wT7sN\/\/\/lIsAdd\/EFq+q1LacC5Y6LeNve8+f8H75hBm\/\/rn+nCC79Fl\/7vk+jEU1voxKYWOon\/f2pLmRDWPS8hCkfdIScYDJF\/FRps2zqUdB+33PnidBcvXkyPPfYYdXV1UUtLS\/T34cOHRzuG1MulMFRu64l8wwcv3OJazzdPPZu2\/OjH9JGXPkGPD\/0\/9MEJf1wjrVsh\/vJ376bDf38znTOnK+pdy8\/Rfd+gN3\/ZE\/3TjzedRo90nU9T5q2mM4b8\/C1RHzqoJty8GO7EU0dEt5KFJOYhCke+r2621ELkHyfcYoorG1H\/scUUnZyTuPu4daOR8lD5oUOH+k0zyr8LYccctx+\/Wwu3uhhNzbaPg+yz9LhZuB944IHoC3PYsGH0i38aS1\/aP5sOnflRWvWJsfTqc28S97jbvj4qmq\/m5\/Wnvh\/1ps\/97MNRL1s8Jx47QCzaR\/f9A5103kcjIWZBXv+JPVGQ\/7F2BB1c91l689ltdOJpA+i8N3fQgGGDIiHvee3s6C7d504YH8VrGXSI\/u+xD9GIQc\/To8c+TBPOejJKg3vsT+x6B51x9s\/prJ+MoH3nvh79++lb99HLz59LL53+Jp35yok06v2n06HHfhr9xh8MF095lZ565OIozLCWHXSg5xI6c9TA6LeXD53b5\/+mVZE\/OORHpGUS\/4033qDnDz5PZ59zNg0YMMAkCsLkSCBE\/k8++0P6y0XXEh+4xO+7eMTIWY54vSQ1+BM30eBPfr6fbd193Lr5bLlt3759O6lrleTfyybc9y3\/ZxrR0rc9inPCsz3P0oIFC2jDD77Spx54cZqj0aCFmw9Kue2v\/4q2vrklwjd9+nS6fuK3ad+vLqaFT86kx559jdb\/4XB69Nqj9KGvnErnvuckosM9dHzDAqILJ9IJl9\/QR7RPe+Yz9OxrQ+i1c6+lfW\/+Nu3\/xRu0+b9foZG9J9M53xhAoz49IBJaTvd3T32e6Olt0f+f+8UbNOLIDvq9AU\/QY6eMpwlnPRUJOg1uoQEDn6fDT7wjsnPg2Uvo1deG01PdH4r+PvyiH9MZQ3qjP798aGgkwOIZ\/s4fUc\/236QzR79Apw16llre85bw\/0fXDVG8p7ovi8Q6Evazfx4Juoivirlp3XKNZ5o+woFAFgLPnfRj+tagz\/YTbu5xv\/722pQs6fuOe\/LbU3VqPnT3cVdNuNt++VU64\/iv2z8TX0x\/rLnawq0uRkiD4mOxgsviNN5rfeBfHqbBy86kESNGRD1ucb3nKb+zjFZ97zn6\/kOH6H1\/P5Balr2D3viN4\/QnP7yRfvbIv9O5t\/037X3hCO194bWop3xJz8einuPVOxbSf+4\/nab93rAapg\/8xmA6fPvL9NLeY\/Sxzw2l0e9\/q+fe73mhh6i5hY7\/7BE69Mwv6aXHn6E9P3iFjgwYFQnzmSNPpvP\/8Ky3euSXHI6E981X34h678cHDqMTjh2I\/m\/6cC9ffc589UTT6JnCHTnyGj3\/\/PM07LxhdDJ63JlYukQOkf8jj2yjjr+e2U+4Xfg1Spy4+7irNlS+8Z5\/MRLhAwcO0CPbttGqVavon7Z+wyhOGX2de4\/bVyFdtoMJ4f7zH1xW28J15KkvRUPdZ314a1SU7z\/0PG2ZuZ9e+ZtB1Lv7O\/TXT\/0l3fu+Lrq9Z3Qk2DyUveJdX4vCbnp9Hk36nd+O\/vyB3zirHwo+OpWfcVMHx64y55vI+OCWvT94JTr05aV9r9O4tsE06v2n1YbqfTHO026Ic6x58suaVoj8Q1uclnYfN7eZYjic6xMWp2V9q+oX31q465c1O0viC5K3ffEqSJMDWIRwt331RDrjg38eGeQ5at7PzcLNPVtxM9jCH7+Dnvjm39GowYPoifd8uibO0Zz23n+gpotuoJPPnpiaaR6eZwEf3XpaNA8tPyzWLNos1Dz\/zD1sMa+emnCDBQhROMrkohD5hybcJvdxy9vBbrzxRuK7J8RdE3HbwUQ91g3B17uOu\/jUJU69y5Vmz1m45eNFJ0+eTDfddBPdfPPN1Nraqt16lZaRPH63PfJUCPdVX9zbZ1GHfNmIuDCkY8P2aBX58Jv\/o9Y7f\/35R+jVHQvoHa3rrFaEc6\/68fWHidOOeuBtg4n\/7cyRAxN743kwKksaIQpHWdhzPkLkX4UGu0x1qAx5cfGpS5wylFXOg5Nwc8Hnz58f7feTVx\/29PTQjBkzqKOjw5t42wAWwv3RG7ZG27rEw2eWnzLyKmq66NNR73j7l39Ik3\/\/T0hetcnbvd7a9vVsvyNSbfIQatgQhaNMvg6RfxUa7DLVoTLkxcWnLnHKUNZMwq2bN5G3DcjzJk1NTWUrb5\/8COG+bPxMuuDrx2u\/yWeWC+H+xIxbaPjn\/6MWRsyF2\/a2Sw2kjpkLUTjqiDfVVIj8q9Bgpzo2sAAuPnWJUzas1j1ude5YFWof+7hdocrCPeorz9SGwLk3zb3u5inPRD3up+78Fsm9cv79lR8uMJ7Xds1fleOFKBxl8meI\/KvQYJepDpUhLy4+dYlThrLm0uMWR+Gpwi0fo9cIPe5vf2Mr\/dX7r46GysUCNQbEws0Lzh6+dQI9v30X\/fn3L4u4sWgfeWpV9Gf1FrGyObfM+QlROMrkjxD5V6HBLlMdKkNeXHzqEqcMZc0k3BxZFJxPHNuzZ09tSwGvSFyyZEl0EpnuUo+yFZ573Czc0z9yE4159wf6zXPzcaTr2n6TBl18ae0ubh4i53ltFvWQjijN23chCkfeDLOkFyL\/KjTYWXxexbguPnWJUzZ21kPlogDqCm7+d+6F84I1vpe7ER4h3P9rxTN0WvddxMPl4uEV4y93T6Nvzr+Rhv3x5XTF7SOj3vYvuqfVFq41QhnLmscQhaNMvgiRfxUa7DLVoTLkxcWnLnHKUNbMPe6yFcI1P0K4F\/3X\/6A3l06iwZ+8qTZczscgvvjwJHro726hsy+ZRB9ZfkI0r83PO96\/ztUk4r1NIEThKJPzQ+RfhQbbpg6pN3ilXXOclnbaPRVp8Yv43cWnLnGKyHuWNJ173FmMliWuLNwn3vVnUbZ45TiLNu\/Z\/tVrj9LX1\/5tdEHHez68MvpdnKhWljI0aj5CFI4y+SpE\/lVosG3qEIRbT6sK9QDC\/Y2txD3uIS\/9P3ru838Q7dXm58hPvkfnzumiNX\/yGl3U+u906ed+1+hkNJsXK+SwIQpHmfwdIv8qNNg2deiGG26gzZs3R1H4\/ojLL7+cZs6cGa0\/WrNmTfTv48ePp7Vr11Jzc\/Nb7d6RI9HRpyIeH661fPny6N857s6dO6NwYh2TYCrypaZnk1+XsC4+dYnjkrci4zgJd9qlI\/V2nisgucfN19SJ6\/14dTkvSOP\/6+7idrWHeL8mEKJwlMn\/IfKvQoM\/E55KAAAgAElEQVRtU4fietyjRo2KxJgfFmmxQ0h3PbJ8eZM6VC4fxMXrmtRjp23y6hrWxacucVzzV1Q8J+FWz8AtKnNFp6sKt84ehLsYL4QoHMWQdEs1RP5xDTYvOq3Co+5yMRkql7fzcm963bp1kaiLrbyyWA8ZMoTmzZtHnZ2dsQuQdbc0FsnWRYRd4hRZBpe0rYXb5PIOl4z4iJMm3Hx++B3vfZKu3nhBZS\/78MGdbYYoHL5Y6+yGyD+uweYtnuJshjL5yDYvvEWVj2kWj61wi+28ql2xWyhOuPnqUD7qev\/+\/VHU2bNnRxc91eNxEWGXOPUoi40NCPfbc9w8VK4+EG6bqmQXNkThsCNUbOgQ+Sf1uH9VgV73Sae29DlbwkW45eOr1RqoDpXLgi3mvNHjLva9FalbC7eoDGJepD7ZLMYKetzFcDVJNUThMOFSrzAh8q9CT8umftgKNw+Vs\/DKi9Vke6pw8zD77t27+\/SuIdw2HnIPay3cbEp8aa1cubIhTkiLw5Mm3OIu7jmPviu6chNPfgRCFI786GVPKUT+EO4XopXhPIwtTraU57i5VsmL1fjvzEyIOf9dnuPW3VPBJ2diqDz7+5mWgpFw605JS0q4UVeVq2WCcKdVH\/ffQxQOd1r5xwyRf2jCzbWGxVWI6axZsxKFmxekqW292paziPNWMt5eNmXKlD5bx1iwx4wZQ+vXr4\/tteddk1186hIn73xnTc9IuLMaKWt8vvnrn277EXGPWjfHLYS7s3dcWYvQsPkKUTjK5KwQ+VehwS5THSpDXlx86hKnDGWV8wDhThDuxzccJh5Oh3DnX21DFI78KbqnGCL\/KjTY7h6vZkwXn7rEKRs9a+FOO6+2ke7jTutxs3BzGO6R48mXQIjCkS\/BbKmFyL8KDXY2r1cvtotPXeKUjRyEO6HHzaLN4g3hzr\/ahigc+VN0TzFE\/lVosN09Xs2YLj51iVM2ekbCnXbEqVooXrjQ1tZWtrL2y09ajxvCXZwLQxSO4mjapxwi\/yo02PaernYMF5+6xCkbRSPhljOdNlRetgIm5QfC7c9bIQqHP9r9LYfIvwoNdpnqUBny4uJTlzhlKKucB2vhLlsBsuQnTbh5YdqL+47RNRsvzGIGcTUEQhSOMlWEEPlXocEuUx0qQ15cfOoSpwxlhXC\/TQDC7a86higc\/mijx80EqtBg29Qhk5PTbNLLY7TV5mQ19pd66YmaXxefusSx4VSPsOhxJyxOu+\/Kn0U+QI87\/6oI4c6fqU2KIfKvQoNt42MIt55WFeoBhDtFuM8aOZCuuH2kzfuCsAYEQhQOAyx1CxIi\/yo02DYV5IYbbqDNmzdHUXjB8OWXXx6dnMbHnfLpZ\/yoJ6OpC5EnT54cXfPJ\/85x+TxzfsSlIoKpyJeannwRCf\/Gd4HL91zI8cUtZHy3t\/zvSSdxuvjUJY4N93qEhXBDuOtRz\/rZCFE4vICOMRoi\/7gGm28BrMKj3qcQ1+Nm8WQx5kc+m1yEb21tre0Kkoe21aFy5jl\/\/nzq6uqK7ucWx6XyhwGfh67eaSH4i7PM1fi6v2OoXF8zSyvc6pec+MIzecHUChsXJ22Om+\/iHtc2mCYtGGpiFmEsCIQoHBZ4Cg8aIv844eZ2YOvK3sKZF21g0vyhfdoqk6Fy+aIQ7k2rQimLddx93HK5ZKFXLyHhcOL3uXPnRh8N8keC\/DsLP+a442tM7sItDqFnk\/LQh02lTfsSS0tL5CFN7CHcaSSL+z1E4SiOpn3KIfJP6nHz7pFGf3haT+512wr3pk2bogtJ1Ee043HCLQ+Hc1zRo9YtRBNXgeouPOG48lWhEO46CnfWyh9337fJakT1ZhsId1ZvFBc\/ROEojqZ9yiHyr8Lcpo2nXYS7u7s7Gkbnm8LURx0qlwVbtLVyOw3htvGWXdjce9x25vuHFuIr3xnLofilS7vknRdP8PzNddddR3PmzKG0+8LTetxLhz4eLUzj4XI8+RIIUTjyJZgttRD5Q7iT7+PmofK0Nla9j3v37t3RfLZ4MFSe7b00jV064eavuIULF9KKFSuiBQ\/iUYfPkwqoLoqICwvhNq0m+YcLUTjyp+ieYoj8IdzJws21SV6spnaY+O+qcMs9dPnubxZztROGxWnu76saE8J924\/oo\/94enQf97Bhw2p8Xn3uTeLFaW1fH0WjW0\/PjzhSigiEKBxlcn2I\/B999FFqb2+nLVu2RO97CI8sprp5ZXUBmTrdqG7FEuuHeHvZlClTIqEXW854bnvMmDG0fv16Wrt2LTU3N9dWlu\/fv594a9kZZ5wR\/Sd66XHbwdg3Ii8HDx6srVxXfSbif+c730n16YEDB6LoPT09DV8PnIRbda4KM2nfXdrL4qPHveHUv4iyNX36dLr22mujP7Nwb7ryVfrQV06lc99zUlq28bslgaNHjxK\/kPyxNGDAAMvYCJ6VQIj8eSj4M5\/5TFDCnbWelD2+EO7777+\/T8dLl+977rmH7r333tpPjfwB5yTc\/NXFX1C8iIFXIop5jjjRtXF+PYWbr+zk88gv+souGjFiROR40evmVaZf+8DuqMc94vcH2RQBYQ0I8MKZ3t7e6CsZwm0ALOcgIfLnRp7XvzRyg51zNWj45IRw\/+u\/\/qtRj5t73du2baNVq1Y1dD2wFm7R2546dWq0SV9dsi8v50+rFfLWMQ7LQy1x2wTSFqfJtkznuIVwT3+suZ\/T+VAGHiq\/euMFGCpPc6TD7yEO1TpgKixKiPxDm+MurPKUKGEXn7rEKVGRo6xkFm4WyaVLl9Jtt91Wm9OQ\/25b4CzbwYStPIR7T\/crdP+VT9OcR9\/VZ2+kbXkQXk8gROEoU10IkX8VGuwy1aEy5MXFpy5xylBWOQ\/Wwq0ei6fb25dFuDlzAqx8Hq58tF4aRAh3GiH\/v4coHP6p\/zoHIfKX25VQFqeVqc4VkReXhWZBCjfD5+FweeUgD3nzakIeOtcdc+fisLQjT5MOZIFwuxCvb5wQhaO+hJOthcifG\/kFCxZEc5x4qkNgwoQJ0VGtpk+wws2A5AVq8s0xWVaUm4LPK1zSHLf4rbN3XF7mkI5EIEThKFMFCJU\/izf\/5\/th\/i+++CKdffbZWJyZ0Rk8emIzghK0cGdkXYroEG5\/bghVOPwR72sZ\/P16Avz98Ydw+2Ofi2UIdy4YnRJBw+WELbdI4J8bSqeEwN8JWy6RghRudTGaSlKd\/86FdEGJJAk3H4fKv\/Oqcjz5E0DDlT9TmxTB34ZW\/mHBP3+mpilCuKWzxAU0CLdp9Qk7HBouv\/4Hf\/D3S8Cf9WCEW2wBE2fSpiHnc2x5hXnZH\/S4\/XkIwuGPPVsGf\/D3S8Cf9WCEW0acNlTuzx32lpOEm49C5WNPr9l4oX3CiJFKAMKRiqjQAOBfKN7UxME\/FVFhAYIU7sJoekgYwu0B+tsm0XD5Y48et1\/24O+XP4TbL\/\/M1iHcmRE6JwDhdkaXS0TwzwWjcyLg74wuc8SghVu+2pPvWb3pppvo5ptvptbW1oaY32bvJwn3fVf+LKogGCrP\/J5oE0DDVQxX01TB35RUMeHAvxiuJqkGK9xccHF2+Pbt26m7uzu64pNPJJoxYwZ1dHQ0hHinCfdZIwfSFbePNKkLCGNJAA2XJbCcg4N\/zkAtkwN\/S2A5Bg9SuNVLRtSzyfM6qzxHP8UmBeGuB2W9DTRc\/tizZfAHf78E\/FkPUrjV+7h1wi1fQOLPPemWIdzpjIoKAeEoiqxZuuBvxqmoUOBfFNn0dIMUbvW+bFW45ctHmpqa0il6DJF05\/Yd732SxrUNpkkLhnrMYXVNo+Hy61vwB3+\/BPxZD1K4Gbd8r+2ePXtqc9ybNm2iJUuWkLhH259rzCxDuM04FREKwlEEVfM0wd+cVREhwb8IqmZpBivcjEdeVS5wDR8+nLq6umis5ihUM6T1DQXhri9v2RoaLn\/s2TL4g79fAv6sBy3c\/rDnZzlJuJcOfTxaUc7D5XjyJwDhyJ+pTYrgb0Mr\/7Dgnz9T0xQh3KakShoOwu3PMWi4\/LFHj9sve\/D3yz9Y4U67dGT8+PG0du1aam5u9uuhFOsQbn\/ugXD7Yw\/h8Mse\/P3yD1a4eeU4F74RxDmpisQJ90v7jhGvKr964wU0uvV0v7WsotYh3H4dC\/7g75eAP+tBCre6j9sf\/uyWIdzZGbqmAOFwJZdPPPDPh6NrKuDvSi57PAh3A9y5jR539opeRApouIqgap4m+JuzKiIk+BdB1SzNIIVbPYDFDFU5Q6HH7c8vaLj8scccq1\/24O+Xf5DCzch37doVXSaycuVKmjhxol8vZLAeJ9zi3zHHnQFuSlQId3FsTVIGfxNKxYUB\/+LYpqUcjHDrDltJgtPoq8qTVpunVQr8bkYADZcZp6JCgX9RZM3SBX8zTkWECka4i4BXhjTTetxzHn0XnTlyYBmyWrk8oOHy61LwB3+\/BPxZh3AXyF7AFSZMzj\/nC0\/4rHTxzJ49mxYtWhSbSwh3gQ7EULk\/uAaWIdwGkAoMAv4Fwk1JOkjh5mHzefPmUWdnp\/ZMchbPrNd6Mtj58+fXzj1X\/67zC9tdvXp1LY4Y3uc5+DjxjhNucd1nZ+84f7Wr4pbRcPl1MPiDv18C\/qxDuDWXiWQV7rhV63zoCz86ERZxWltbqU3aopYm+E\/859P04J++QuqQOIS7+JcKwlE84yQL4A\/+fgn4sx6McKcdcaq6YNmyZX0E1MZFoqfMAi2vWGfYLN42p7Xx6veFCxfSihUrtKMDQrjV1eMQbhuPuYWFcLhxyysW+OdF0i0d8HfjlkesYIRbhpU2VJ4VbJzYpvWe44bPk4bt44R76629xOLNPXE8xRBAw1UMV9NUwd+UVDHhwL8YriapBincJmCyhMlLuE3muIVwf3D1KTSq9TQaNmxYlPX\/\/NtDtP3+gxDuLI5MiYuGq0C4BkmDvwGkAoOAf4FwY5I+cOBA9EtPTw+1t7fTli1bqKWlpf4ZycHiCcePHz+eQzq5JZGHcIuh\/b179yYOrQvh\/tagz9JzJ\/2Ypk+fTtdeey396K6j9My33qApG0\/LrVxIqC+Bo0eP0sGDB6OPpQEDBgBPnQmAf52BK+bAv\/7877nnHrr33ntrhiHcOfogq3CbijZnWe5xnzjqcCQi\/B\/3uH\/yDy\/TzO+PybFkSEomwH7q7e2Nvngh3PWvG+Bff+ao\/36Zc4+b\/9u2bRutWrUKPW5Xd\/BiszVr1tSi877rWbNm0cyZM6PV47aL02xEWxZudXHaQ3P30Yv7jtE1Gy90LRripRDAUKHfKgL+4O+XgD\/rmOMugL3LdjDOhq1oQ7gLcJ5FkhAOC1gFBAX\/AqBaJAn+FrByDgrhzhmoSE6AFaelmawo5947h7PZLha3qhw97oIcKyWLhqt4xkkWwB\/8\/RLwZx3CXSD7tCNP5QNZxG1l+\/fv1+Yo7rhUCHeBDsRQuT+4BpYh3AaQCgwC\/gXCTUkawq0BJM9bDx8+vHYEqT83xVuOE+77rvwZnTVyIF1x+8gyZrsSeULD5deN4A\/+fgn4sw7h9sc+F8sQ7lwwOiUC4XDCllsk8M8NpVNC4O+ELZdIEO5cMPpLBMLtjz0aLn\/s2TL4g79fAv6sQ7j9sc\/FMjvwOx89ldTtYBgqzwVvYiIQjuIZJ1kAf\/D3S8Cf9SCFO+2s8qy3g9XTnUK4eS57XNvgmuk73vtk9PdJC4bWMztB2YJw+HU3+IO\/XwL+rEO4C7jWs57uhHDXk3ZfWxAOf+wxVO6XPfj75R+McNfzWs96uhTCXU\/aEG5\/tPtbxoeTX2+Avz\/+wQi3jDhtqNyfO+wtxwn30qGP06T5QzFUbo\/UOAYaLmNUhQQE\/0KwGicK\/saocg8YpHDnTtFjgknCrc57e8xmJU2j4fLrVvAHf78E\/FkPWrh1p5WV\/cAVtapAuP29PBAOf+wxx+qXPfj75R+scIuCL1u2jNra2mpe4BXlS5YsobgjRv26q791CLc\/j0C4\/bGHcPhlD\/5++Qcp3HG3dwlX8JGnfGb48uXLqampya+HUqzrhPulfceIt4Ope7tLXZAGzByE26\/TwB\/8\/RLwZz1I4ebFaXxf9tSpU\/v0toUbGn0fN4S7Pi8UhKM+nOOsgD\/4+yXgz3qQwo0et78KVyXLEA6\/3gR\/8PdLwJ\/1IIWbcVd5jhs97vq8UBCO+nBGj9svZ\/AvH\/9ghZtdUdVV5Xu6X6H7r3wac9wFv28Q7oIBpyQP\/uDvl4A\/60ELtz\/s+VlmBz74p6\/QtFW\/WzurXAj3nEffRWeOHJifMaTUhwCEw2+FAH\/w90vAn3UItz\/2uVgWwv2Reb9TOyUNwp0L2tREIBypiAoNAP6F4k1NHPxTERUWAMJdGNr6JAzhrg9nnRU0XP7Ys2XwB3+\/BPxZD0a4xRawnTt3GtEeP348rV27lpqbm43C+wqUJNydveN8ZSsIuxAOv24Gf\/D3S8Cf9WCEW0UsFqZ1dHT0Ozlt9erV1NXVRWM1V376c5Xesk64H99wmB6au48g3MV6C8JRLN+01ME\/jVCxv4N\/sXyTUg9SuKu2j5sXp8lz3BDu+rxQaLjqwznOCviDv18C\/qwHKdxVOzkNwu3nBYJw+OEurII\/+Psl4M96kMJd9R731lt7iXvdvB0MT3EEIBzFsTVJGfxNKBUXBvyLY5uWcpDCzVCSTk5r9DluCHdatc\/ndzRc+XB0TQX8XcnlEw\/88+Hokkqwws2wdCvNJ0+e3BC3ggln6xanQbhdXgX7OGi47JnlGQP886Rpnxb42zPLK0bQwp0XxLh0BFzxu8kd32oc9b5w1RaEu2gvxqePhssfe7YM\/uDvl4A\/6xDugtgz2Pnz59e2lal\/15kVW9RWrlxJEydOrJ2lrm5Zk+MK4X7\/Jy+iK24fGf3EW8Fe3HeMrtl4YUGlQ7IQDv91AMLt1wfg749\/kMLNQ+Q7duygyy67rBDycYvfbrnllsjeokWLtHb59\/379\/cZque7wbu7u2OH79mB93\/safqjj0+CcBfiTfS464zV2ByEwxhVIQHBvxCsRokGK9wzZ86MAOlOR2Phvf3222nWrFlOJ6eJuXMWaO45i4dhszjbnMimE3O1xw3hNqrruQdCw5U7UqsEwd8KV+6BwT93pMYJBincTEeeS5bnkVko16xZQ1mOPOUh74ULF9KKFSv6nL5mMlwue04dOtd5Vdfjvu\/Kn0VBMVRu\/B44BUTD5YQtt0jgnxtKp4TA3wlbLpGCFW6mJ4a0N2\/e3Afm7NmzY4ezTahnFW55tXvaKvc44T5r5MDa0LlJnhHGngAaLntmecYA\/zxp2qcF\/vbM8ooRtHDrxNtk5Xca\/KzCLafPc9zr16+PHV4Xwn3x+86nK1aNpGHDhtH6T+yh0887EcKd5qiMv6PhyggwY3TwzwgwY3TwzwjQIfqBAweiWD09PdTe3k5btmyhlpYWh5T8Rznh+PHjx12yIYbFOe5dd91FDz74IHHvO8swOaeVp3CnHc8qhPuVE35O3zvldpo+fTq1bPsEnXbeiTTxxkEuWBDHkMDRo0fp4MGD0cfSgAEDDGMhWF4EwD8vkm7pgL8btyyx7rnnHrr33ntrSQQl3PJQtDosLoYgTMVbFn+myenxojZe\/JbH4jRT4eYe9yULTohE5LsdR2l06+n0vr8akqWOIG4KAZ5q6e3tjb54Idz1ry7gX3\/mskXwrz9\/7nHzf9u2baNVq1aF1eNmMbzzzjtp7ty51NTU1I9+1lXlLtvB4uKkLVDTzXHf8d4naVzbYJq0YGj9a1ZAFjFU6NfZ4A\/+fgn4sx78HHdR6AVYMWdusqJcjSPEnPO4fPly7UcGhLsoD6anC+FIZ1RkCPAvkm562uCfzqioEBDuosgqW87YjLrwTXcgi3rkadoKdw6\/5qpH6KqrrqotRkOPu0CnSkmj4aoP5zgr4A\/+fgn4sw7h9sc+F8tCuD\/0+5Nr+7aXDn2cJs0fiqHyXAjHJwLhKBhwSvLgD\/5+CfizDuH2xz4Xy3HCzeeW8zw3nuIIQDiKY2uSMvibUCouDPgXxzYtZQh3GqGS\/w7h9ucgNFz+2LNl8Ad\/vwT8WYdw+2Ofi2UIdy4YnRKBcDhhyy0S+OeG0ikh8HfClkskCLcGo7w3e\/jw4bWrOXMhnnMiqnC\/tO8Y8eK0qzdeEO3lxlMcATRcxbE1SRn8TSgVFwb8i2ObljKEO41QyX+HcPtzEBouf+wxVO6XPfj75Q\/h9ss\/s3UId2aEzglAuJ3R5RIR\/HPB6JwI+DujyxwRwp0Zod8EINz++KPh8scePT6\/7MHfL\/8ghZuPPN2xYwdddtllfunnYJ0duPLjm+gPhk2nOY++i\/Z0v0L3X\/l09OczRw7MwQKSiCMA4fZbN8Af\/P0S8Gc9WOHmS0D4Wbt2LTU3N\/fxQNazyuvpTgh3PWn3tQXh8McePT6\/7MHfL\/8ghZuRy0eLLlu2jNra2iJPiBXlpreD+XXfW+VAj9uPFyDcfrgLq+AP\/n4J+LMerHAzcnGJB9\/BLT9p54P7c1d\/yxBuf96AcPhjjx6fX\/bg75d\/0MKtE2\/1IhC\/7km3rgr34xsO00Nz91Fn77j0yAiRiQCEOxO+zJHBPzPCTAmAfyZ8mSIHLdzyQSt33XUXPfjgg8S970YZJmfPQ7gz1f9MkdFwZcKXOTL4Z0aYKQHwz4QvU+QghZtXlfPitJ07d5I6LC6ANIp4Q7gz1f9MkdFwZcKXOTL4Z0aYKQHwz4QvU+RghfvOO++kuXPnUlNTUz+Ajbaq\/KZP3kWXHp0bDY9vvbWXeLict4PhKZYAGq5i+aalDv5phIr9HfyL5ZuUepDC7Q93\/pbZgRDu\/LmapIiGy4RScWHAvzi2JimDvwmlYsIEI9xieHzRokU0ceJEY5oMiOfCdfu9jRMpMCCEu0C4KUmj4fLHni2DP\/j7JeDPOoQ7hT2E21\/lLLtlCIdfD4E\/+Psl4M96cMLNC9JsnzIvVFN73LwV7MV9x+iajRfaFhPhLQlAOCyB5Rwc\/HMGapkc+FsCyzF4MMKdI7NSJQXh9ucONFz+2GOo3C978PfLH8Ltl39m6xDuzAidE4BwO6PLJSL454LRORHwd0aXOSKEOzNCvwmwA+e1LaX\/+doXoy1gm+fuo7NGDqQrbh\/pN2MBWEfD5dfJ4A\/+fgn4sw7h9sc+F8sQ7lwwOiUC4XDCllsk8M8NpVNC4O+ELZdIEO5cMPpLBMLtjz0aLn\/s2TL4g79fAv6sQ7j9sc\/Fsirc9135NI1rG0yTFgzNJX0kEk8AwuG3doA\/+Psl4M86hNsf+1wsQ7hzweiUCITDCVtukcA\/N5ROCYG\/E7ZcIgUr3OIu7uHDhxOfprZhwwZasmRJBHXy5Mm0fPly7TnmNtQFXBHH9spQjj9\/\/nzq6uqisWPHak1DuG08km9YNFz58rRNDfxtieUbHvzz5WmTWrDCzceYcuH5KNNDhw7RjBkzIsHmi0cWL15MQtBtYMphVdE1EWE5vjii9eDBg6nC\/RfTPk1Tf\/lVunrjBXT\/lU\/TpPlDMVTu6jiLeGi4LGAVEBT8C4BqkST4W8DKOWiQwi16262trdTW1hYJeHt7O4keMfe+169f73w+udqbFz7jjwV+uIef9oi7wvkDIq3HrQo3bwXjeW48xRJAw1Us37TUwT+NULG\/g3+xfJNSD1K4RW926tSpkXCzUK9evbomkFmFO+5CE9Nzz0W4WbNm0dKlSyHc\/t6PRMtouPw6BvzB3y8Bf9Yh3G1t0e1fYti8ubk5EvLu7m7nee5du3bRwoULacWKFX3mpk2Gy+WPitGjRxvNcYseN\/e0+axyHjIf3Xq6v1oViGUIh19Hgz\/4+yXgz3qQws24Waz3798fCez1118fXfXJQ9hqb9zFNVmEW\/5o4AtRTBan8TD\/p179x2hue+vKXmr7+igIt4vjLONAOCyB5Rwc\/HMGapkc+FsCyyH4gQMHolR6enqi6d0tW7ZQS0tLDinXP4kTjh8\/ftzWrJiH3rx5M4nbv5qamnJZmOYq3Go8kx66+PJi4d5+8np6z+tT6UNfOZXOfc9JtkgQ3pLA0aNHiRcPDhs2jAYMGGAZG8GzEgD\/rASzxQf\/bPxcYt9zzz1077331qIGJ9wu0EzjuAi3umCObdkKd\/MHXqUXvn8aXff9MdF55XiKJcA+6+3tjb54IdzFstalDv71Zy5bBP\/68+ceN\/+3bds2WrVqVXg97ryQi9XfIr3Zs2cTLyqbOXNmNPTOQ\/DiSVqcxmLPW9J4+F73LFu2LFpIpz5yj5tXkj++4XB02ciZEO68XBybDoYKC0ecaAD8wd8vAX\/Wg53jZuRiPpvnknkP90033UQ333wziW1irm7JYzuYbY+77ZdfpYvfdz7t7X6VOnvHuWYd8SwIQDgsYBUQFPwLgGqRJPhbwMo5aLDCLQ9Db9++vbaKnCf9uefb0dGh7eGa8lf3hpsMe+t606aL01i4ecj2pX3HINymTsoYDg1XRoAZo4N\/RoAZo4N\/RoAZogcp3Op8srr9K+t2MOGPtCNP0w5kMRF7YUMIN9vmoXI8xRNAw1U84yQL4A\/+fgn4sx6kcOsOYJH3bWc9gKWe7pSF+4zj50Zz2xDu+ngAwlEfznFWwB\/8\/RLwZz1I4dZdMCILt9jjncdFI0W7Vjhw9pC76Vd7z4JwFw1cSh\/CUUfYGlPgD\/5+CfizHqRwM255DnrPnj21Oe5NmzZFt4TZ3uTly4WqcI9qPY2u2Xihr+wEZRfC4dfd4A\/+fgn4sx6scDNyeVW5cEHapR7+XKW3LBy48N33R3u4Idz18xCEo36sdZbAH\/z9EvBnPWjh9oc9P8uqcPNebj6zHE\/xBCAcxTNOsgD+4O+XgD\/rEG5\/7HOxLBz4N1O+Ts+sOxl3cedC1SwRCIcZp6JCgX9RZM3SBX8zTkWEClK4eYh83rx51NnZ2ef2LgG4EVeVb77\/3+7BpTMAAAumSURBVOiHK47TpAVDccFIEW+KJk00XHUCHWMG\/MHfLwF\/1iHcY8f2o9+Iwt3Ih837q\/7ZLEM4svHLGhv8sxLMFh\/8s\/HLEjsY4ZZvAzMBFnc2uEnceoapggPryStPW2i48qRpnxb42zPLMwb450nTLq0qtPvW13qmDZXbIfQbugoO9EvQ3ToaLnd2ecQE\/zwouqcB\/u7sssasQrtvLdxZoZUpfhUcWCaeNnlBw2VDK\/+w4J8\/U5sUwd+GVr5hq9DuQ7jb2xv6XtZ8q3T9UkPDVT\/WOkvgD\/5+CfizHqxw6w5fkd0wfvx4Wrt2LTU3N\/vzjoHlKjjQoJilDALh8OsW8Ad\/vwT8Wa9Cu+\/U45bPI+djTnfv3k2LFi2iXbt20cKFC2nFihXarWL+XKW3XAUHlo2paX4gHKakigkH\/sVwNU0V\/E1J5R+uCu2+tXCrt4MxhHXr1pG4VIS3gwkhzx95vilWwYH5Eqlfami46scaQ+V+WYN\/ufhXod3PLNzcy166dCnddttt0dC4+vdyuaxvbqrgwDLzTcobhNuv58Af\/P0S8Ge9Cu2+tXCLPd2tra3U1tYWXTYin6QG4fZXIRvJMoTDr7fAH\/z9EvBnPUjhZtzq6Wg85z1mzJhIyPk3+X5uf+5Jt1wFB6aXspwhIBx+\/QL+4O+XgD\/rVWj3rXvcAre8QI174TNnzqSdO3dSo6wo53JUwYH+qn82yxCObPyyxgb\/rASzxQf\/bPyyxK5Cu+8s3FnAlSVuFRxYFpa2+UDDZUss3\/Dgny9P29TA35ZYfuGr0O5DuHEAS35vhEVKaLgsYBUQFPwLgGqRJPhbwMo5aJDCnXZWOW4Hy7mWVTQ5NFx+HQv+4O+XgD\/rEG5c6+mv9jW4ZQiHXweCP\/j7JeDPejDCjWs9\/VWyqlqGcPj1LPiDv18C\/qwHI9wy4rShcn\/usLdcBQfal7ocMSAcfv0A\/uDvl4A\/61Vo97E4DYvTvLxBEA4v2GtGwR\/8\/RLwZz044eZT0WbMmEGTJ0+OLhXhR0DgPw8fPpy6urpyuWBETpfTfuCBB2jixImJ3ua95WvWrOkThvMqzlFXI1fBgf6qfzbLEI5s\/LLGBv+sBLPFB\/9s\/LLErkK7b9zjFqLd0dERnZAmi7YQVQHERGSTwHM68+fPr30EqH\/XxVWPYjVxbBUcaFLOMoZBw+XXK+AP\/n4J+LNehXbfWLjlk9Kampoi6vxvDEG+e5v\/jR\/RI7d1jxBg7r3LaaSlK24t4zhpPXORpyo40JZvWcJDOPx6AvzB3y8Bf9ar0O4bCbeuNxsnsFn3cccJMMNm8ZY\/EmTXu1xuUgUH+qv+2Szz1a933303fepTn6KWlpZsiSG2NQHwt0aWawTwzxWnVWJVaPeNhFu9g5sp6f6N\/z2rcLMAL1y4kFasWNFnrjxtuJztLlmypI8Dk+a3OWAVHGhVY0sUGOz9OgP8wd8vAX\/Wq1D3nYVbzHmvXLmyz9C0bkjdxkWuws12N2\/eXJsXFyMCbDttcRrPyaPXZ+Ol7GF7enqovb09WnQI9tl52qYA\/rbE8g0P\/vnytElNsN+yZUvDtj1Gws1Q1Dlm7uGuXr26zypy0QvnOWbXOW5X4dY5Lu7jQoRlBy5YsIC2bdtm43eEBQEQAAEQaGACEyZMoHXr1jVsCYyFWx6q5h7S4sWLo0KL3qzo4e7duzd2HtqEUp7CHTecL+eDxZv\/wwMCIAACIBAGAdawRh7pMxZudqe8t1qePxb\/Pnv2bKuetrrvmuPPmjUruttbXR2etjhNV91MhDuMaopSggAIgAAIVIWAlXDXo9Au28Hi4sT13utRDtgAARAAARAAgSIIlE645Z69fLCLfCCLDoR6+Ese8+1FAEeaIAACIAACIJCFQCmFWx2W57+rp7HpDmRRj0m1HbrPAhJxQQAEQAAEQKAeBEor3PUoPGyAAAiAAAiAQKMRgHA3mseQXxAAARAAgaAJBCncYjEbH9jCT9oJa0HXkAIKrzvlbvz48Zm2ERaQzUolKer8tGnT+p3l73ITX6Xg1KEwcfzFWpydO3f2ycWyZctqlznVIXuVNaG2Nbrp00as\/8EJt7oCPW5FemVrcgkKlnZhTAmyWLksiK2X6loRl5v4KgenDgWK44+dL8XBVw8J0y1YbtT6H5xw6848x8tT3Mujpuxy\/Wr9clc9S2qPThZul62X1SNUbImS+LNll\/Mpis1xNVKPa2d0B4nZ3kRZBkLBCTd\/hXV3d\/c5vxxiUr+qyA3ZvHnzqLOzs88lMvXLQTiWhGiMGjWKrrvuOpozZw7Jdwu43sQXDsFsJU3jz6lze8Q3hbkeEZ0th+HFljtpQ4YMye2wr3qTDE64dcO0GC6vX7VT55PYMua3i+evO7c\/z+OFiy9BY1vQ8VfX2ogSYn67OF\/Lt1ceOnTI6SbK4nJnnjKEm4gg3OYVJmtIsVhEHrLljykW9Li71rPaRHwiCLffWqDjL\/fIxZ0PIlxHRwcWp+XsMnWOu5E\/XCHcEO6cXw\/75HCmvD0z2xgQblti+YZPu6lQtib3Cpubm\/PNSKCp6S7BgnA3UGXAUHn5nIURj+J9AuEunnGSBRvh1i2g9Zv7xrYed3MlhLuB\/IrFaeVzFoS7eJ8kDdXmcRNf8SVobAsQbj\/+S7puupEXZwY3VI7tYH5eIGGVRzz279\/fZ1V\/3AvkN6fVsp60OKoRt8M0mnfiRjxmzJjRZ6U\/l0vXuWi08pYhv0mizflr5O2QwQm3cBY7jheE8LN48WJSG68yVLwq5kFdfKP6o6mpqYrF9l6muB6feqsehmmLcVUcf3VhpuqPYnITRqomi14btf4HJ9zylxaOPPXzAotGjHve\/ODI2eL9kDRU24hHPhZPLF8LSfzFqWrConq6Xb45CSM1tY1RSy0zbsT6H6Rwh1F1UUoQAAEQAIEqEoBwV9GrKBMIgAAIgEBlCUC4K+taFAwEQAAEQKCKBCDcVfQqygQCIAACIFBZAhDuyroWBQMBEAABEKgiAQh3Fb2KMoEACIAACFSWAIS7sq5FwUAABEAABKpIAMJdRa+iTCAAAiAAApUlAOGurGtRMBAAARAAgSoSgHBX0asoU0MQEGe079y5MzG\/y5Yto9GjR9P8+fOpq6uLxo4d67V8tmfLc\/h58+ZRZ2en97x7BQfjIJATAQh3TiCRDAhkJdAo95LrrsZNKzsfK8nx1q5dS7hjOo0WfgeBZAIQbtQQECgJgUYQblcBFpfJtLa2UltbW0mIIxsg0JgEINyN6TfkuoIEkoRbvbVLXI96ySWX0Be+8IWIBt9wx0Pphw4dovb29hoh3aUVfHXkkiVLamFmz55NfC930hN3DaLuQgedTVfRr6CrUSQQyEQAwp0JHyKDQH4EbIV7zZo1JARXiCrfeCfftsYCvXr16j5z4yz6HE7Mlwu7EydOTBTvpDulOzo6aj3puKtBG2FEIT9vIiUQKI4AhLs4tkgZBKwI2Ao3C6Q8Z6wTaVVs1fvQRQZN7oHWCbLN\/d1xPXYrSAgMAiBAEG5UAhAoCQFb4eb7zJcvX05NTU1RCUyEWxeG45r0hjnu+vXr+3wsyCvjTYbbxRC\/nO+S4Ec2QKBhCEC4G8ZVyGjVCdRLuOW5bZUpbz2LWzzGoqv28kV8\/o2H7sUTJ+IQ7qrXYpSvHgQg3PWgDBsgYECgXsKt9poNslbr0ZvEFQvfdOIN4TaljXAgEE8Awo3aAQIlIVAP4dYtMOPix\/27jMZ0PlvMZXNceUgcc9wlqWjIRsMTgHA3vAtRgKoQqIdwMyt1VXmc0KpcdeKuW9QWtwDOZB69Kr5EOUCgSAIQ7iLpIm0QsCBQL+HmLLns4xaiz\/+X93wL8ZaLqpsrxz5ui8qAoCCQQADCjeoBAiBgTMB0uFxNECenGSNGQBBIJQDhTkWEACAAAjIBnFWO+gACfglAuP3yh3UQaDgCuB2s4VyGDFeMAIS7Yg5FcUAABEAABKpN4P8Dl4JFZ9xZhB8AAAAASUVORK5CYII=","height":238,"width":395}}
%---
%[output:9791a366]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"ans","value":"            Peak: 0.4323\n    SettlingTime: 4.7690\n       Overshoot: NaN\n"}}
%---
%[output:15857bb6]
%   data: {"dataType":"textualVariable","outputData":{"name":"xb_max","value":"0.0100"}}
%---
%[output:6adab004]
%   data: {"dataType":"textualVariable","outputData":{"name":"theta_max","value":"0.0873"}}
%---
%[output:0e6468c6]
%   data: {"dataType":"textualVariable","outputData":{"name":"disp_percent","value":"10.6016"}}
%---
%[output:195c8ce3]
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
%[output:8be810fa]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAe4AAAEqCAYAAADTbrP0AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnV+IHdmdmM8YEbY3JoQeLwipZ1Ag7SRP2pgEa0X8GPYlmiGG0Jp+mFm5PQhsxg9xW1KbsY1fpFZb8yCGsNFD0zPzoG2Rh8BoISRvIy+NRcIwHQIh1kNkIcYa5NGLycoKy3ZySj53qktVt+rrfzq36mswHt37+9Wt+n51fl+dulV1X9ja2toK\/klAAhKQgAQkMBEEXlDcE1EnV1ICEpCABCRQEFDc7ggSkIAEJCCBCSKguCeoWK6qBCQgAQlIQHG7D0hAAhKQgAQmiIDinqBiuaqTQ+Dx48fhwoUL4ebNm+H48eNhdXU1TE9Pjzbgxo0bYWlpqfa9Llt5+fLlcO3atXD27Nlw\/vz5Lim7jknr3LSgU6dOheXl5RC3fWFhIWxubjZ+5vXr18OJEyfCnTt3wpkzZ8Knn346dv0uXboU5ubmdr0NLkACfSCguPtQRbchOwJlcceVqwpWcSvu7HZaV2hiCCjuiSmVKzpJBKrijuueZpnxv3cr7ufBous6P3r0aDTjpjPl9Blp9j41NfU8NtXPlEDWBBR31uVx5SaVQJ24yzIaJ8HqKem6U+11p8rLwkzcmgRY\/YwuglXck7o3ut59I6C4+1ZRtycLAtXvuB8+fFh8j5sEWSfBOtmnjTly5EhYW1sLs7OzxUtVcddJO+VWT9On3Cqotu\/LFXcWu5YrIQHv43YfkMB+EChLOArx2LFjxcVoScAff\/zxMxen\/eIXvwjz8\/PF6qTT6uXllGfPVXGni7xibhJ8Wl5Z+um18iw+ST8eXJQPDqpcxl2cVv6McQcRcZnjDhA8Vb4fe6PL7BsBZ9x9q6jbkwWBqrjffPPN0fe+ZZGXBdp0pXidgNtm3NUZeoLSdmX4uFPmijuLXcuVkIAzbvcBCewHgaq44y1b5Rn1N77xjfDzn\/98dDtYvAgr3T5WlWfdbLpO8k23VtUdHDRtcxdx133nXl6eF6ftxx7lMiXwBQFn3O4NEtgHAnXijh9T\/X55r2bc1U2ofl9e\/W59J1dt+x33PuwoLlICOyCguHcAzRQJtBFoEnd1VlwWd\/lUNP2Ou+50etusv\/oQlLhNXb7jdsbdVn3fl8D+ElDc+8vXpQ+UQJO4I46yoMsS3M1V5V1zx8V5VflAd1Y3e+IIKO6JK5krPAkExom76+NQ03bWndZuupCteiq+6SK1alybtMsHHM64J2EPdB37TCBbcZcv5IkFKD91qqkg1ateuzSjPhfXbZOABCQggf4RyFLcUdqLi4vb7kct\/7uuDFHa77777ignXdkaf8jgoH6EoX+7h1skAQlIQAK5EchO3Ok0YjzFVxZuPLUX\/+oknHJOnjy57ReEqgcAucF3fSQgAQlIQAKUQHbiTjPlKOg4W05\/UcJR3tWfRxy3wfEK3nPnzoWVlZXRoyIpIOMlIAEJSEACORHITtxNst3J7DmePl9fX0eyz6k4rosEJCABCUigSqC34vY7bnd2CUhAAhLoI4Feijt9533v3r3W2fb9+\/f7WFe3SQISkIAEGgjMzMxMNJveiZtK+wc\/+EG4ffv2RBfRlZeABCQgge4Evv71r4ef\/exnYVIFnp24d3NxGpF2LHG6VzwW8OjRo92rPtDIeIBz9erVYoeXV7edQGbdOKUoeTFeMVpmjFnidevWLcXN0DVH7+R2sLg0Ku2yuCe5gHvFvcty0oGOvLrQehojs+6s5MVYpWj3McatD7yym3GXB296WlqXK8rjrWIxjtwu1ocCsl12d9F3794N7733Xnj77bfDoUOHdrewgWTLjBVaXoxXjJYZY9aHvp+luMvyTiWpPvK0\/ECWpt8hbsr1SJXt6Cn6d7\/7Xfj1r38dXnrpJcXdEaHMOoL6fZi8GK8YLTPGTHEzXtlF96GABwnVBsFpy4wxkxfjpbg5rz70\/Wxn3LwcPKMPBeRbvfMMmypnJzPGTF6Ml+LmvPrQ9xX3\/HzwYqtuO79NtRuncpTMGDN5MV6Km\/NS3JxZVhl9KOBBArWpctoyY8zkxXgpbs6rD33fGbcz7s57vk21M6pRoMwYM3kxXoqb81LcnFlWGX0o4EECtaly2jJjzOTFeCluzqsPfd8ZtzPuznu+TbUzKmfcHFWR4T7GwcmMMVPcjFd20X0o4EFCtUFw2jJjzOTFeHmww3n1oe8743bG3XnPt6l2RuWMm6Nyxi2zHRJgaYqb8couug8FPEioipvTlhljJi\/Gyxk359WHvu+M2xl35z3fptoZlTNujsoZt8x2SIClKW7GK7voPhTwIKEqbk5bZoyZvBgvZ9ycVx\/6vjNuZ9yd93ybamdUzrg5KmfcMtshAZamuBmv7KL7UMCDhKq4OW2ZMWbyYryccXNefej7zridcXfe822qnVE54+aonHHLbIcEWJriZryyi+5DAQ8SquLmtGXGmMmL8XLGzXn1oe8743bG3XnPt6l2RuWMm6Nyxi2zHRJgaYqb8couug8FPEioipvTlhljJi\/Gyxk359WHvu+M2xl35z3fptoZlTNujsoZt8x2SIClKW7GK7voPhTwIKEqbk5bZoyZvBgvZ9ycVx\/6vjNuZ9yd93ybamdUzrg5KmfcMtshAZamuBmv7KL7UMCDhKq4OW2ZMWbyYryccXNefej7zridcXfe822qnVE54+aonHHLbIcEWJriZryyi+5DAQ8SquLmtGXGmMmL8XLGzXn1oe8Pesb93b\/4n+E\/\/KePwscX\/1WYmZnhe8DAMmyqvOAyY8zkxXgpbs5LcXNmWWUoblYOmyrjZVOVFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GK7toxc1KYoNgvBS3vDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GK7toxc1KYoNgvBS3vDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GaxT96NGjsLCwEDY3N4vXzp49G86fP995aZcvXw7Hjh0Lc3Nz23Li69euXdv22qlTp8Ly8nKYmpp6ZvmKuzPyItAGwXjJTF6cAM9wXDJmipvxKqKTtE+fPl2It\/rvtkXeuHEjLC0thUuXLm0T9+PHj8OFCxfCyZMnnxF60zIVdxvt7e\/bIBgvxS0vToBnOC4ZM8XNeBXRUbzr6+thdXU1TE9PF69FkHG2XH6tuugk5ps3bxZvVcWdDgDizP3EiROd1kxxd8I0CrJBMF6KW16cAM9wXDJmipvxKqKjoONf+dR4m3STtO\/duxeuXLkSFhcXQ5qxp1W4c+dOuHjxYnjnnXdGBwRtq6e42wg542aEno22qTKC8mK8PDjkvBQ3ZNZ0OpucLm+KTafQy6s07vvtGKe4WQFtqoyXTVVenADPcFwyZoqb8Qr7Ke44k4+n0dfW1sLs7Ozos+Iqtl2c9uF3\/jjMzMyEw4cPwy0aVrgNgtdbZoyZvBgvDw6783rw4EERfP\/+\/TA\/Px9u3bpV9P1J\/Htha2tr66BWfD\/FXbcN8fT5mTNnitPrdd97pxn3l\/9qpUh\/\/fXXwxtvvHFQOCbuc548eRIePnxYHOAcOnRo4tb\/eaywzBh1eTFeMVpm3Zi9\/\/774YMPPhgFK+5u3PZ1xl23Cm2n4JO4\/\/zVPwpHjx4thOSsu7mY8cDrs88+K45SFXe3nV5m3TilKHkxXjFaZt2YxRl3\/N\/t27fD1atXnXF3w\/Y0aicXp5WX3yZjEut33KRy3sfNaD2N9tQvoyYvxst9jPPyO27ObMe3g6WPqhN3OgV\/5MiRbVerx1Pl586dCysrK8X33tU\/xc0KaFNlvGyq8uIEeIbjkjFT3IxXEZ3EG79zjreEkRl0Ob96O1gqxvXr14vvs6ufU7eqipsV0AbBeClueXECPMNxyZgpbsZrFN32yNN4a9fGxkbt1eDjRJ8Kkj6o7VGqipsV0AbBeClueXECPMNxyZgpbsYru2jFzUpig2C8FLe8OAGe4bhkzBQ345VdtOJmJbFBMF6KW16cAM9wXDJmipvxyi5acbOS2CAYL8UtL06AZzguGTPFzXhlF624WUlsEIyX4pYXJ8AzHJeMmeJmvLKLVtysJDYIxktxy4sT4BmOS8ZMcTNe2UUrblYSGwTjpbjlxQnwDMclY6a4Ga\/sohU3K4kNgvFS3PLiBHiG45IxU9yMV3bRipuVxAbBeClueXECPMNxyZgpbsYru2jFzUpig2C8FLe8OAGe4bhkzBQ345VdtOJmJbFBMF6KW16cAM9wXDJmipvxyi5acbOS2CAYL8UtL06AZzguGTPFzXhlF624WUlsEIyX4pYXJ8AzHJeMmeJmvLKLVtysJDYIxktxy4sT4BmOS8ZMcTNe2UUrblYSGwTjpbjlxQnwDMclY6a4Ga\/sohU3K4kNgvFS3PLiBHiG45IxU9yMV3bRipuVxAbBeClueXECPMNxyZgpbsYru2jFzUpig2C8FLe8OAGe4bhkzBQ345VdtOJmJbFBMF6KW16cAM9wXDJmipvxyi5acbOS2CAYL8UtL06AZzguGTPFzXhlF624WUlsEIyX4pYXJ8AzHJeMmeJmvLKLVtysJDYIxktxy4sT4BmOS8ZMcTNe2UUrblYSGwTjpbjlxQnwDMclY6a4Ga\/sohU3K4kNgvFS3PLiBHiG45IxU9yMV3bRipuVxAbBeClueXECPMNxyZgpbsYru2jFzUpig2C8FLe8OAGe4bhkzBQ345VdtOJmJbFBMF6KW16cAM9wXDJmipvxyi5acbOS2CAYL8UtL06AZzguGbNBiPvx48fhwoUL4bXXXgsnTpxghDKPVtysQDYIxktxy4sT4BmOS8ZsEOJ+9OhRWFhYCKdPnw5zc3OMUObRipsVyAbBeClueXECPMNxyZgNQtxpxn3y5EnFzfaP3kXbIHhJZcaYyYvx8uCQ8xqEuCOWO3fuhHPnzoWVlZUwOzvLSWWa4YybFcamynjZVOXFCfAMxyVjNghxp1Plm5ubjXSOHz8eVldXw\/T0NCP4nKMVNyuADYLxUtzy4gR4huOSMRuEuBmSyYpW3KxeNgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjEblLjL33WfOnUq\/OQnPwk\/\/elPwyRfba642Q5vg2C8FLe8OAGe4bhkzAYj7rihi4uLYW1tLXz88cdhY2MjLC8vh\/v374czZ86Et956ayJvFVPcbIe3QTBeiltenADPcFwyZoMQd\/U+7hs3bozEPTU1Far\/Zgifb7TiZvxtEIyX4pYXJ8AzHJeM2SDEXX1yWp2419fXvR2M7TsTGW2D4GWTGWMmL8bLg0POaxDiTjPuI0eOhPPnzz8zw758+XL49NNPi1PncQY+SX\/OuFm1bKqMl01VXpwAz3BcMmaDEHdEkjb0+vXr4Ve\/+tXoVPmHH34YlpaWQnx9En+ARHGzHd4GwXgpbnlxAjzDccmYDUbcEUvdE9TiLDxesDapj0FV3GyHt0EwXopbXpwAz3BcMmaDEjdDMxnRipvVyQbBeClueXECPMNxyZgpbsYru2jFzUpig2C8FLe8OAGe4bhkzBQ345VdtOJmJbFBMF6KW16cAM9wXDJmgxB3+bvtg\/wVsAQ3lYRcAJeuhH\/ttdfGXjSnuNkOb4NgvBS3vDgBnuG4ZMwGIe6IJN7yde3atUY6ey308pPa4oVv1X+3lSmtb5vsFXcbye3v2yAYL8UtL06AZzguGbPBiLsJS5qNx\/f36ve4q\/eNp8+OMo5\/8V7ytvVJvx2uuNkO3RZtg2gj9Oz7MmPM5MV4eXDIeQ1S3Hfu3CmeTx4fuhL\/9vqWsHQwEAVdvjc8wo7ybjpASHkvv\/xy+Na3vhW++93vhitXrniqnO\/XjRk2VQ5TZoyZvBgvxc15DUbc+y3rMvr4WefOnQsrKyvb7g8np8vT+ipuvlOPy7Cpcp4yY8zkxXgpbs5rEOLej9Ph41A\/D3F\/+J0\/DjMzM+Hw4cN8LxhQhk2VF1tmjJm8GC\/F3Z3XgwcPiuD4q5bz8\/Ph1q1bRd+fxL8Xtra2ttpWPH3vfPPmzSJ0ry9Ge94z7i\/\/1UqxCq+\/\/np444032nAM9v0nT56Ehw8fFgc4hw4dGiwHsuEyI7RCkBfjFaNl1o3Z+++\/Hz744INRcO\/FXcay3xJ\/HjPuP3\/1j8LRo0cLITnrbh4EsfafffZZcZSquLs1C5l145Si5MV4xWiZdWMWZ9zxf7dv3w5Xr17t\/4y7Cct+nEbf6cVp1Vl7vIDO77i77dBdozyN2ZXUF3EyY8zkxXjFaJkxZoP6jjvdYlWHaC9Pne\/mdrC0bl6cxnbkrtE2iK6kFDcn9TTDfYyTkxljNjhx7\/WtX024yz8jGm8JI1eUx2UqbrYjd422QXQlpbg5KcUts50SYHmDEDdDsnfRbY88HfdAFsW9d3UoL0lxc64yY8zkxXh5loLzGpS4q\/dyR1wHNQPnpemWER95euOj\/xE+eftPJva2gG5bujdRNlXOUWaMmbwYL8XNeQ1G3GlDL126FObm5kakbty4EZaWlkLbo0U52oPJUNyMs02V8bKpyosT4BmOS8ZsEOJuulgsoYqnrOPjT5eXl8PU1BQj+JyjFTcrgA2C8VLc8uIEeIbjkjEbhLjT7VmnT5\/eNttOqOKse319fc9+ZISVYHfRipvxs0EwXopbXpwAz3BcMmaDELczbrZT9DnaBsGrKzPGTF6MlweHnNcgxB2x+B033zn6mGFT5VWVGWMmL8ZLcXNegxF3RONV5XwH6VuGTZVXVGaMmbwYL8XNeQ1K3BxP\/hl+x81qZFNlvGyq8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMei\/u9ICVhGVSH7TSVFbFzXZ4GwTjpbjlxQnwDMclY9ZrcUdpv\/vuu2FtbS3Mzs6OLk576623au\/nZujyiFbcrA42CMZLccuLE+AZjkvGrLfiTvdunzx58plHnG5sbEzkU9LqSqu42Q5vg2C8FLe8OAGe4bhkzHor7qanpcUNjo84XV1dDdPT04xWhtGKmxXFBsF4KW55cQI8w3HJmA1S3IuLi6PT5wxXftGKm9XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLlpxs5LYIBgvxS0vToBnOC4Zs96Le3NzsxOR48ePT+QFa4q7U3lHQTYIxktxy4sT4BmOS8ast+JmGCY3WnGz2tkgGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GK7toxc1KYoNgvBS3vDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GK7toxc1KYoNgvBS3vDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFK25WEhsE46W45cUJ8AzHJWOmuBmv7KIVNyuJDYLxUtzy4gR4huOSMVPcjFd20YqblcQGwXgpbnlxAjzDccmYKW7GK7toxc1KYoNgvBS3vDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZrxG0Y8ePQoLCwthc3OzeO3s2bPh\/PnzY5f2+PHjcOHChXDz5s0i7tSpU2F5eTlMTU2N8i5fvhyuXbu2bTl1cSlAcbMC2iAYL8UtL06AZzguGTPFzXgV0Unap0+fDnNzc8\/8u26RSdpHjhwpBF\/9d8xJr508ebJYbpc\/xd2F0hcxNgjGS3HLixPgGY5LxkxxM15F9I0bN8L6+npYXV0N09PTxWsRZJwtl18rLzq+v7i4GNbW1sLs7Gzx1p07d8K5c+fCyspK8Vo6IIhiP3HiRKc1U9ydMI2CbBCMl+KWFyfAMxyXjJniZryK6Cjo+Fc+Nd4m3Sj7jY2NbafGqzPsKPKLFy+Gd955Z3RA0LZ6iruN0Pb3bRCMl+KWFyfAMxyXjJniZrwaT2dXT59XF1sn++rp8ij3paWlbanjvt+OgYqbFdAGwXgpbnlxAjzDccmYKW7Ga1\/FHeUeL1xLp9OT2OMqVi9iS6utuFkBbRCMl+KWFyfAMxyXjJniZrz2Vdx1qxJPn585cyZcuXKl9nvvJO6\/\/LOXwszMTDh8+DDcomGF2yB4vWXGmMmL8fLgsDuvBw8eFMH3798P8\/Pz4datW0Xfn8S\/F7a2trYOasWbrvzei1PlddvQttwk7r\/3X57eivb666+HN95446BwTNznPHnyJDx8+LA4wDl06NDErf\/zWGGZMeryYrxitMy6MXv\/\/ffDBx98MApW3N24FVH7dXHabsR97V9+KRw9erQQkrPu5mLGA6\/PPvusOEpV3N12epl145Si5MV4xWiZdWMWZ9zxf7dv3w5Xr151xt0N29Oo\/bgdLIokPpwl3eed1qd6y1h1Pf2Om1QuBE9jMl4xWmaMmbwYL\/cxzsvvuDmz0f3W8V7reEtY2+nsdEQZxRz\/4oVm8a8q6lSM69evF99np+Wmz6lbVcXNCmhTZbxsqvLiBHiG45IxU9yM1yi67ZGn4+7bHvfI01SQ9EFtj1JV3KyANgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLlpxs5LYIBgvxS0vToBnOC4ZM8XNeGUXrbhZSWwQjJfilhcnwDMcl4yZ4ma8sotW3KwkNgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLlpxs5LYIBgvxS0vToBnOC4ZM8XNeGUXrbhZSWwQjJfilhcnwDMcl4yZ4ma8sotW3KwkNgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLlpxs5LYIBgvxS0vToBnOC4ZM8XNeGUXrbhZSWwQjJfilhcnwDMcl4yZ4ma8sotW3KwkNgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLlpxs5LYIBgvxS0vToBnOC4ZM8XNeGUXrbhZSWwQjJfilhcnwDMcl4yZ4ma8sotW3KwkNgjGS3HLixPgGY5LxkxxM17ZRStuVhIbBOOluOXFCfAMxyVjprgZr+yiFTcriQ2C8VLc8uIEeIbjkjFT3IxXdtGKm5XEBsF4KW55cQI8w3HJmCluxiu7aMXNSmKDYLwUt7w4AZ7huGTMFDfjlV204mYlsUEwXopbXpwAz3BcMmaKm\/HKLjqK+y\/+64Pw3\/\/tbJiZmclu\/XJbIRsEr4jMGDN5MV4eHHJeipszyypDcbNy2FQZL5uqvDgBnuG4ZMwUN+OVXbTiZiWxQTBeiltenADPcFwyZoqb8couWnGzktggGC\/FLS9OgGc4Lhkzxc14ZRetuFlJbBCMl+KWFyfAMxyXjJniZryyi1bcrCQ2CMZLccuLE+AZjkvGTHEzXtlFJ3F\/8vafhJen\/yC79ctthWwQvCIyY8zkxXh5cMh5KW7OLKsMxc3KYVNlvGyq8uKTL9eqAAAMXUlEQVQEeIbjkjFT3IxXdtFJ3B9+55+Gf\/EP\/35265fbCtkgeEVkxpjJi\/Hy4JDzUtyc2b5mPHr0KCwsLITNzc3ic86ePRvOnz\/f+JlJ3Of\/9Fg4\/6f\/YF\/XrQ8Lt6nyKsqMMZMX46W4OS\/FzZntW0aS9unTp8Pc3Fyo\/rvugxU3K4dNlfGyqcqLE+AZjkvGTHEzXvsafePGjbC+vh5WV1fD9PR08VmxQJcvX972Wnklkrhf++eHw7977Z\/s6\/r1YeE2CF5FmTFm8mK8PDjkvBQ3Z7ZvGVHQ8a98ajzNuuNrJ06ceOazFTcrx927d8N7770Xvv3tb\/ts947oZNYR1O\/D5MV4xWiZMWaKm\/Hat+jHjx+HCxcuhJMnTxanydNf2+lyxc1K0ocdnm3x7qNlxhjKi\/FKZxbn5+fDrVu3PKDugK8P+9gLW1tbWx22NeuQ3Yo73sP9l3\/2UtbbmMPK3b9\/P8QGcf36dRtEx4LIrCOo34fJi\/GK0TJjzBKvST7QGbS4L\/\/n\/x1+9h\/\/W\/jbP\/xK+NJf\/6ao\/pf++vPRXnDoN\/+r+O+\/+co\/Kt6PcX\/7hy+G+Hr675ST3k\/JMa68rPLr8b\/TMp7+\/4uj5cWc+O\/yX1pOeR3S55Zj47KerusXy0j\/XV1u\/Hd5eekz0nqn5Zbzxw2P6jqn2Oq6xNfLn5HWI\/7\/\/3355DMfUbfedXyaPr+6Hk3Lq25ndR2rdY3syqyrnFMtyhtU91p5HylvQ10N036T9snytiWudftNlXl536nb3+pqV11u235Sty3VcVO3zHHr9nfubdTuI+UxGGtS\/vfTbf9K8VrXdWoag9UxXe4FTctP61LdscvjK+WWx295fxs37pp4te3nTcusGwflzxhXn7rP7LoddcxTvZv2m7rxVObdNJ7SOv3m3\/+bNrTZvj9occeqxKOv+L9fPflyIbK\/+co\/DvcePQ4vT08Vr8e\/+Fvd8bW6hhbj0nvxv+Nfyi\/npLi03DbRlJdb3nuaXid7WBwkcZvi9nVdj\/LnVgfZbtepzI1sR13suKbbtuxxEmvL3S2DanMsf15ar\/JnVJtWl3Xfi3Vs45Deb\/qspgPJury92C+q7Kr7cVexVLc7LaftwKeuf7RJsyvjGFftP3W9p+vyuuwf1e099gf\/J9z93d8d+xE7qWPddqW+lXpsdfu7jIFyD\/\/Bv\/5nXdFkF9cLcUeqO7k4LbtquEISkIAEJCCBFgK9EfdObgdz75CABCQgAQlMGoHeiDtdQR5v+4q3f7VdUT5phXJ9JSABCUhAApFAb8QdN4Y+8tRdQAISkIAEJDBpBHol7kmD7\/pKQAISkIAEKAHFTYkZLwEJSEACEniOBAYp7vTAlps3bxboT506FZaXl8PU1NPbuYb0F6\/GP3bs2LYnznX92iE9gSjxig9mqT5aNl40uLS0VIQcOXIkrK2thdnZ2YlCXN3Ouv2ly9c0Q+EVi3vnzp1w5syZ8Omnnxa1rvulPpnVD4PE7sqVK9vGk7y28yr3lvTO8ePHt\/02RV+ZDU7cSdpRIvEituq\/J8oou1zZtONfunRp7KNi6y70ixJaXFwcibj677hq1Sv966783+Um7Ht6km06KEn7S\/zgdLDX5ZfphsKrfNBX\/aW+dOHouJiUE2OGxCztyOVJRflA2H3s2aFedwtwOarPzAYn7jrBxCPcc+fOhZWVlYmbDe7EXNUzDlVxt91aF89MxGfDp4OftA7lgVQn+6ZH0+5kGw4qp645VGdE8np2JtT2S30yq9+Dy7PIsrjltZ1Xl17SZ2aDE3cs5sbGxrZT4112goMSxX5\/TtrWe\/fuhXgqLs6ay7Oc+PltD7P56le\/GhYWFoozFuVT4+WfUf38889rD4bq+O\/3Nu\/18qu\/OievdsLVJiqzZ5mlCcSbb74ZLl68WIzPNL7ktZ1XHIPf\/\/73ww9\/+MPGyVafmQ1O3HXFHOrpcjIrLsd+7Wtfq5Vy+WxGFHfdb6FP4unyaostb2d8dGzbL9MNnVd1P+vyo0BDY1buQd\/85jeL6wOSuOX17EFO9XqRGFH+frvvzBT3\/\/8Nb8V9evQdd993+Pa54fiI6nfcMVpx1zMrfyUzpKa6k32serZKcY+nmL5SKH+dECcKkePq6mpxoXGfx6XiVtzbTpUr7vENIzaHeDdCujpeXt005VmKZk7Vr16q11C4j3Xbx8pndl555RXF3Q3bZER5qvyLOnmqnO2zVWnHbJtqN4bls1rf+973et1UuxH5IqrakxQ3Jfg0fkj72OBm3EO\/OK08JJqe597nizp21hKeXrBXnmmXlyOvdqrVr6Nk9pRZ9T7jKsl0\/7u83MfKBAYnbm8HGz\/jju\/2+TaK9uH\/bMQ4acurnld88Er5oUbVg0T3seY9se4BLPLaziuOyaZ9LN3t0mdmgxN308VF1XuSd9LgJy2nacbd5ZfWqg8maXoAy7vvvjv6PngSryiP61zehroay2s7lSSet956a3TRY\/nCoenp6dFMc9yv+Q1lH6vuU3Xidh8bv4+NezBSH\/exwYm7\/F3I0B95Ou6nT\/v6qEBycNV2GrP84Bp51TfW9MhTHxPbfc\/zkafdWFUfqzukfWyQ4u62WxglAQlIQAISyI+A4s6vJq6RBCQgAQlIoJGA4nbnkIAEJCABCUwQAcU9QcVyVSUgAQlIQAKK231AAhKQgAQkMEEEFPcEFctVlYAEJCABCShu9wEJSEACEpDABBFQ3BNULFdVAhKQgAQkoLjdByQgAQlIQAITREBxT1CxXNWDIRAfz3nt2rXGDyv\/BvDBrFH\/PqX8+9PxEajVR6JWt7jt\/TpC8XG1d+\/eDfHZ1f5JoE8EFHefqum27AmBnUhiTz54IAup\/v503Ow25m3v16Fr+snVgWB2M3tMQHH3uLhu2s4I7EQSO\/ukYWbV\/URlG\/O295tIVmf2wyTuVveNgOLuW0Xdnl0T6CqJNHN85ZVXwurqavEzg+n3k9NsL\/2QTVypulPs1R9K+NGPfhQ++eSTkH6trumHYOp+tazLZ6ZfaIunjxcXF4t1jn\/lH0xJAGPs0tLSiGeKqf5SVQoY96M1KabpBzTamJffn5qaChcuXCh+H73ur8zZWfeuh4MLyJCA4s6wKK7S8yXQJpGqqB4+fDj66dL4XhLYyy+\/PPpN6upPVMa49FoSYlm86QCgq7i7fmaScfmXlNJrZeFVf4O8+lOddYzqftq1Wsmmn3ZtY972ft32lw9A1tfXi4Or+H26fxKYdAKKe9Ir6PrvOYFxF6eV5dY082z6De+qfOK\/44x3eXk5xFlk\/EuCjGKNs+Ku4u76mXVx1c9omhXH3I2NjWJ9Nzc3w\/z8\/LazCHXbUy1OU0zbBYFxOcePH6+VbzrguXfvXu37XQ4o9nwncoES2EcCinsf4broySTQNrurzrhPnz4d5ubmipeTROJ\/l4WcZtjx9PTa2lp48cUXw8LCQijnlvPJqfKZmZni1HHbZ87Ozoa6GW9V3F1El7azbT3Le0A1p\/xeG\/Nx7zcdtKTlNx2ITObe6VpLIATF7V4ggQqBNol0EXfT969RdPsl7rbPbBP3iRMnill+mwjT9pcPAn75y18W35nHbYufU\/fXdIYixrYxb3q\/7iuI6mcrbod43wgo7r5V1O3ZNYE2iXQRd5qJNq1Mk8S6zmTLck0z7rbPjOsybsZNxV0W4kcfffTMaf\/qtu\/1jLv6tUITa8W96yHhAjIjoLgzK4ir8\/wJ7EbcTXKse30333FXLx5ruuir+noXcTeJrnoKvfy1QPx+Oc7Wo\/zH\/Y37jjsuv+kCsmpNxl2MVv38Lqf+n\/9e5xpIoDsBxd2dlZEDIbBbcddJpXpVdkRZnTHWXVVedxo5nR5Op93jqemun9lF3Okz46n3dOq76QxBWpemC8equ8xeXVXetUbjDqQGsju7mT0koLh7WFQ3aXcEukph3H3LXe6pjmuZlhGv0o5\/8Urt3\/72t6P7uONr1WXFK85fffXV8OMf\/3jbd8pdPrOruJPw6u7jLtMd9711XRXqnppWd3BSzS3X5PPPPw9nzpwZ3YNejS3fk+593LsbC2bnSUBx51kX12qgBMZ9D5wjkp18f9zltrG92lafnLZXJF1OTgQUd07VcF0GT2DSxF2+tzvdi95WxKZZd1sefd\/ZNiVm\/KQQUNyTUinXcxAEJkXc6WlrXb\/brhbvIGbC\/jrYIIbMIDdScQ+y7G60BCQgAQlMKoH\/Byq\/SiXpFZs2AAAAAElFTkSuQmCC","height":238,"width":395}}
%---
%[output:22ddab3e]
%   data: {"dataType":"textualVariable","outputData":{"name":"tau","value":"0.0114"}}
%---
%[output:0f6679e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"num","value":"1"}}
%---
%[output:97a12f32]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"den","rows":1,"type":"double","value":[["0.0114","1.0000"]]}}
%---
