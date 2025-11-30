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
% M1 = [(I_b/(r^2)) + m_b, m_b*l;
%       m_b*l,            I_r + m_b*(l^2)];
% 
% % M2 = [[1,  1+(g*m_b)],
% %       [1-(m_b*g), 1 - (M_r*d*g)-(m_b*g*l)]]
% M2 = [1,             (g*m_b);
%       (m_b*g), 1 + (M_r*d*g) + (m_b*g*l)];
% 
% % M3 = [[m_b*l],
% %       [-m_b]]
% M3 = [-m_b;
%       -m_b*l];
% 
% % --- 5. Calculate A and B ---
% % This is the numerically preferred method for A = inv(M1) * M2
% A1 = M1 \ M2;
% 
% % This is the numerically preferred method for B = inv(M1) * M3
% B1 = M1 \ M3;
%%

A = [0 1 0 0;... %[output:group:229a734b] %[output:78eafac0]
    (m_b*g*r^2)/(I_b*l) 0 (m_b*g*r^2)/(I_b) 0 ; ... %[output:78eafac0]
    0 0 0 1;... %[output:78eafac0]
   20 0 (M_r*d*g)/I_r 0 ] %[output:group:229a734b] %[output:78eafac0]
B = -[0 ; 0; 0 ; ((I_w/R) + (M_w+M_r)*R)/(I_r)] %[output:5675b116]
C = [0 0 1 0] %[output:0708fda7]
% C = eye(4)
D = [0] %[output:46cdf202]
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
R = 1e-4;
%%
Q = [1000 0 0 0 ; 
    0 100 0 0;
    0 0 1000 0;
    0 0 0 100];
%%
[K,S,P] = lqr(A, B, Q, R);
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

figure; %[output:6309dbbc]
plot(t, y); %[output:6309dbbc]
xlabel('Time (s)'); %[output:6309dbbc]
ylabel('States [x, xdot, theta, thetadot]'); %[output:6309dbbc]
legend('x', 'xdot', 'theta', 'thetadot'); %[output:6309dbbc]
title('Closed-Loop Response with LQR (ODE45 + Safety Constraints)'); %[output:6309dbbc]
grid on; %[output:6309dbbc]

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
Qdiag = diag([1000,1,1000,1]);   % larger on states you want faster
Q = Qdiag;
R = 0.005*eye(1);     % measurement noise covariance
N = zeros(4,1);               % no correlation
[L,P,E] = lqe(A, G, C, Q, R, N);   % or simply lqe(A,G,C,Q,R)
%steady‑state Kalman gain L, the error covariance P, and a vector E
%containing the eigenvalues of A-LC
Ke = L %[output:19a15f7d]
%%
x0_o = [0.02;0;deg2rad(7.5);0] %[output:973f237d]
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
Q_min = 0.5 * eye(nb);     % nb × nb  process noise covariance
R_min = 0.005 * eye(na);   % na × na  measurement noise covariance
N_min = zeros(nb, na);     % nb × na  cross-covariance (typically zero)

% Design reduced-order Kalman filter
[L_min, P_min, E_min] = lqe(A_bb, G_min, A_ab, Q_min, R_min, N_min);

Ke_min = L_min %[output:63b81036]
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
x0_omin = [-0.002;0;0] %[output:3c41c227]
%%
%[text] ## FFT
% y_noise  = out.y_noise.signals.values;
% ti = out.y_noise.time;
% save('noise_data.mat','y_noise','ti');
%%
data = load('noise_data.mat') %[output:10205beb]
y_noise = data.y_noise;
ti      = data.ti;
Ts = ti(2) - ti(1)
Fs = 1/Ts; 
N = length(y_noise);
Y = fft(y_noise - mean(y_noise));        
f = (0:N-1)*(Fs/N);          
P = abs(Y).^2 / N;           
figure; 
plot(f(1:N/2), P(1:N/2));
xlabel('Frequency (Hz)');
ylabel('Power');
grid on;
title('Noise FFT');

%%
%[text] ## Low pass filter
fc = 8;
tau = 1/(2*pi*fc)
num = [1]
den = [tau   1]

% 10 Hz is far above the 1.5 Hz interference
% the robot's physical angle dynamics are active around 3–6 Hz → fully preserved
% it significantly reduces high-frequency measurement noise
% it barely introduces delay (a few milliseconds)
%%
%[text] ## Butterworth filter
fc  = 8;        % cutoff
Fs  = 1/Ts;     % sample rate from your model
Wn  = fc/(Fs/2);
n   = 2;        % 2nd order

[b,a] = butter(n, Wn, 'low');


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:78eafac0]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["68.7533","0","24.5250","0"],["0","0","0","1.0000"],["20.0000","0","39.4511","0"]]}}
%---
%[output:5675b116]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["0"],["0"],["-12.4258"]]}}
%---
%[output:0708fda7]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":1,"type":"double","value":[["0","0","1","0"]]}}
%---
%[output:46cdf202]
%   data: {"dataType":"textualVariable","outputData":{"name":"D","value":"0"}}
%---
%[output:63c095d8]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":1,"type":"double","value":[["0","0","-12.4258","0.0000","854.3165"]]}}
%---
%[output:09e92a9e]
%   data: {"dataType":"matrix","outputData":{"columns":5,"exponent":"3","name":"den","rows":1,"type":"double","value":[["0.0010","0.0000","-0.1082","-0.0000","2.2219"]]}}
%---
%[output:7f3d3531]
%   data: {"dataType":"text","outputData":{"text":"\nG1 =\n \n            -12.43 s^2 + 4.415e-14 s + 854.3\n  ----------------------------------------------------\n  s^4 + 8.882e-16 s^3 - 108.2 s^2 - 1.137e-13 s + 2222\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 45 49 50 46 52 50 53 56 32 52 46 52 49 52 53 101 45 49 52 32 56 53 52 46 51 49 54 53 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 56 46 56 56 49 56 101 45 49 54 32 45 49 48 56 46 50 48 52 52 32 45 49 46 49 51 54 57 101 45 49 51 32 50 46 50 50 49 57 101 43 48 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:6335ea38]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:47b75a85]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:1b06f103]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0","rows":4,"type":"double","value":[["0.0100"],["0"],["0.0873"],["0"]]}}
%---
%[output:84041628]
%   data: {"dataType":"text","outputData":{"text":"Simulation completed without violating limits.\n","truncated":false}}
%---
%[output:6309dbbc]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAfsAAAExCAYAAAB7+AsAAAAAAXNSR0IArs4c6QAAIABJREFUeF7tfQ2QXcV15kGWMQMCpBFY8jADAjyydzdBDqaK2bFVUEmokCqkDa4lMxpcqMayM+VFizfR30gYMC6jkYRUiQQO0cJkkGLrJ3aiGJGkstryBjlWUALEQ5KqGAVbQpNBQoAgAguBjLbOFf3S03N\/uu\/tvve8N9+rotC813369Pf17e+e\/j3nzJkzZwgfIAAEgAAQAAJAoGEROAdi37DcomJAAAgAASAABCIEIPZoCEAACAABIAAEGhwBiH2DE4zqAQEgAASAABCA2KMNAAEgAASAABBocAQg9g1OMKoHBIAAEAACQABijzYABIAAEAACQKDBEag7sT958iT19\/fT7t27a9TMmzeP1qxZQ01NTdF3O3fupJUrV9LAwAB1dXVVRuHrr79OixYtisofHByk5ubmcb4cOHCAent76dOf\/vSYOlTmtFGwwjLOnzlz5iTWS4r\/Ev2waRc2acy6ueZZu3Ytbd68uWampaWFhoaGqL29vfZdUf45\/0MPPTTOrlk2F7ht2zbq6OiolR2Xhn80\/VT1Hh4ejqU7qZ0q+3o\/oZ7H0dHRVFyKtqu4fsysf1YZuq9x3GXlL\/K7TT9cxH5VeVWbcOXCxV+FHefRdcvFRp60dSX2cQ+iqrTe2CH2eZpCfJ60zp5zmC9a\/kpuXEtJovz000\/TU089RStWrCBX4Wa0bPNkiWNfX1\/kA3+K8K\/KYQFX9lzKThJ71TJUh5xlM07sGeuenp7IlC72+vdxfUtSq3Tpc+KE0qxTVus3bdiKPed74IEHaOHChWNe6rLK03+37YddbBZNyz5t2bKF7r777lrQl8dmXrHXn12bclU7C\/lSYfpRN2KvN269M9Ibnvre5cGzISZvGpvOt14ie3OURNXt2LFj46K2vHhN5Hyqzao2bNN2TLxs86gOzRRBXehUJ5T0LNnwr\/LqHZpL2Ukdr\/pevWhy35A1gqZjZQql3rbz9h0u+dQzf+mll9ZGxlSd9L4t7XlQdXj22Wetnz8bzrKeQZd+OMuWr99Vu60q8DCfXZt6KS4uv\/zy0qL7uhH7NFFksrdv314DLe7BM9\/+496EzUjCbDymjazfH374YXr00Ucj7osO49v4n5VG4XLPPffQj370o9pUSNp0R1Inpj\/0qjM3O1FTTOIimrgoMss\/M7Iwy1E8PvbYY\/S9732vVk+9rCxfmDMzykvCKa7jVT6yHTU0ruzdcccdpIacuV3s2bMnmnbSI8mNGzfSN77xjeir9evX09KlS6M8aRGcjdhndfhmx+XCf5yg6mLkWnaS2Jt4T58+3UnsVZ3mzp1LP\/jBD2qRfVybsI2Y84g9TxWkRXZmf6TaOePMLzf6tIVq20ltNm7047bbbqMf\/vCHY6YQs14iXPph9jOrTyr6rMaNxDCmPELGU1Tcl\/AzprDmtsLTpvo0jf5cm20uyz9z5Eu1F667WY7OdRbONi8MrmnqRuwVqTZvvuaDlzbspAhIGjI0Iy1zXtCMLuLmDdPmtm0iexv\/bdKkDckmCVlSJ2aKWWtr67i1FNwYVd15PYW51kI1VlW2jX9xDzfb0TvltOFf5pt9yvIlqRxbnPS6mFHyV7\/61doLSJbY88gJR4B6u0pqTzZinxUFme3xiSeeiF3\/Evcyo3c+cb7Ylq0iXn5R5g47ThD1Tnn27NnjxE\/3JW4kkJ\/bWbNmjalb0nRAEt5JbUSVnSTkNkPwSc8C1+VLX\/pSrNjfcMMNtakJvf7cZm+66aZxefiF4ZVXXiH9hSyrP3Lph236pKLP6hVXXDGuzrrYx71Am310XN+RpQtsl9McOnRo3Iv617\/+dbr33nvHvFCYfRT\/nXfKwFXkVfq6EXuXoRJToMxhPxYelUaJ9aZNmxI7FgYrrnydLE7Dc4B6x6DyFBV7F\/\/10QazjqrjdvExa87WjCjiyufOZv78+ZHApg07xuGlf8cjJevWrYuEUu\/AzSFQEy\/mRom7jS9xoxZZQqp3gnfdddeYlwn2VX3H9dejdjXikzSMzx2TesHIioyzfNRHK5KGPE0b5qiD2dEkvXzHiUaW2Jtl24h9kpDFib3i9aWXXoqiPVU3ha85vK6\/oMa95OUV+7iIV38xjlvIa2Kn2rR6nvSXbXMtA6fl+qoRAX3qzewrs0YoXPrhPP2W67Oa1KbjyjbbbtxznhTZq+fF9I8XgJuYZL0wKT+ysM4r6kn5Gl7sk0TG7Dhfe+212DdEtTo47Q2UOwL+mDsAsjpO7ig\/97nPpa7GTxru0f1\/5JFHiIesTSE16\/jcc88l+pg0954m9nrkkpbOXEuhv23rq7\/Tpl\/YvziRZFtmJ61EIm0uNmn4jVeipy34ShrW1blWw+5Tp06NqnrRRRfRF77wBbrzzjujIdPly5fT4sWLo9+yxF7nJWvoz4fYm1imiX3aEHScsGeJfRKPaZG9LvY6nnEdnuI8az2CntdW3PJ23OYLg\/kSlvS7KfZqKiNuZFG1WZVGb1O6MN133310\/\/33W72QZ42w2vRb\/Oz\/6Z\/+aRRk5X1Ws8Q+7iUtbSdEktjb+Kcwsd2pYNu2fIl+3Yh93KIWBULanL2t2KvtRkkPl4r844CfCGJvRj8852W7atvsGEyhTVtYqb+wuIp93MuI+fDH+RI3TJr0gqK3B32EQXVgBw8ejEYibrzxxmjYL0mckiJ7XcB8iH3c6AC3efZdj3aV6JjD+HELYuOeiThhzxqZMDHImrNnXBlTNYyfJvZpK+A5H9f3uuuuIxa8uHadJW55xV7vw\/SRQXPUISuyzyv2eptatWoVrV69OnUbsG0\/bL6MqP7VbANK7PM+q7zLI66txbUd82VYH7kxh+2T\/uZ6JY2GmG0ka3slxD7hNSZpFagehSWJhs1wktqjr4pPioqTHnrV4FyGyPWoNG2fvY3\/6mXEdRg\/zu+4yCZpL7L6Pitqi6M1KYrUMdT9cx3Gt+lAlF+6L3o5Lltj9BdFFU3pI0ZmhKWLkzkXGhel+xB7LtNcER+3niLtBUwX\/KQ1DElDmapsc4RExy6to9U7W5fV+DZi\/8UvfpG+\/OUvR02Co864ofGkKMtF7OM6eb39Pfjgg7Rs2TKKG243h5PThvFNX5Paj\/nCm7Zg16Ufdum38j6r\/ILKz5h5Vkmc2KdNlRURe5t1DHGcu7QZH9F93UT2ujDqKynjIi4TxKyFImkLtrIW4JkdeN4FenF14rrxg3fttdeOW9mp6q0aaVYdeToibajdduGZKlfZUvVPWqDH6fXoKw6fpGF+vYGbLxVm409bZGO+jcctVlL2klY1q9\/TtvfoL55xK6fTxMl8UYgbxbAV+6TDZbIWm9o8S0n86wfxcJqkKYW0KRLOF7cWI6mjM+emk+qt2qB+YE9chKa\/COll+t7Slfas6us79IPDzPbHf5trYJLWEMTNN+tY6\/7Y7D5I81\/Pb9MnpYly1noM1VbMcvQFevpLRNoaCx9iz5hy\/X\/7t387elmL++j+YIFexiuMzXxI2tyv6hBstt6ZC+uytpGYv7tsvUsTe14EklW23sEm1VHHhYeX1elpRbfeJXUmZiebxZ2tf+bDnbT1Li1ayPKFfc+aT41rquZiQb0cc6Ed51dz9ia\/f\/AHf0C8FkNP40vsld9ZK6H1F0Szjej1ihPDLF\/jyjZHUZL8y3ouszpZ84VFr5vZLnwLvSo77qUnbtcA9wvcV33zm9+kP\/zDP6zNp6uXa3OdTlab1V\/49brZLGjTcbV5fmz6pCyx574vqyzzd+ZT9W9pbYrxVjsYFPZJc\/ZpfYnJJafN2uKX9XxkyGCun+smstcbqe0wk0Ika74tF3J1mKnsYSNXiKT751qfek2vz+HHrQy3rZe5IM42H9KVj4ASubS+tXyvGrdE9YKAQ3UMjjmS4xXMvO2KP+rf5tChepNcsmQJ8WKTuN8bt\/lk10y6mEr3LxthpNARUB2aflwuEJKHQNaZCfI8rn+P1AiMy5qgorWui8ieRWDfvn21E\/L4LZQPxIi75IYbLq8o3bBhQ+zFM0UBq+f80sVUun\/1zH1VvjOncRfhVOUPyh2LgD5VUqbwTGQe1BA+Y4CLcIyWwA2SP+oyDfNvPbmvIciJ3BhRdyAABIAAEGgsBOoisjcjeY4WeAGGEn+dEnPFuc01rCMjI43FKmoDBIAAEAACXhHgRZH1\/Gk4secXA17BqoZHzL9NsljoeZvE\/v3765lH+A4EgAAQAAIBEbj++uuJz0CoV9GvG7G3HcY3udYX98Ut2FMLJZjEyy67LGBTgek0BPhli\/eWg4fq2wm4qJ4D9gA8yOBB52Lv3r0Q+5C0mMP2aQv04sQ+bcGeEvt6JjEk9mXZBg9lIZ1dDrjIxqiMFOChDJTtymgELuoisrfdeqdWOXZ2dkYr9dXffChF3Pw+09wIJNo1V9mpwIMcfsCFDC7AgwweGkUn6kLsGeykQ3WUoC9YsIB4P6\/rCVh4oGQ8ULzg8vHHHye+633y5MkynJqgXoALGcSDBxk8QOzl8FDIE4h9Ifi8ZX7nnXfo5Zdfpra2Noi9N1TzGQIX+XDznQs8+EY0v71G0Im6iezz05SesxFIDIVNmXbRsZWJdnpZ4EIGF+BBBg+I7OXwUMgTiH0h+LxlRsfmDcrChsBFYQi9GAAPXmD0YqQRdAKR\/dNPU09PD2E1vpdnIrcRdGy5ofOeEVx4hzSXQfCQC7YgmSD2QWAt12gjkFguYmFKQ8cWBtc8VsFFHtTc82Sd3Hn69OloHcvHPvYxrGNxh9cqh+0BOY2gE4jsEdlbPRShE0FgQiNsbx9c2GOVNyVO7syLnN98tqfiQez94l6JtUYgsRLgPBcKgfEMaAFz4KIAeJZZcXKnJVABk6kTCm2mcBtBJxDZI7IP+DjZm4bA2GMVOiW4CI0wDvMKj3B2CS4C7pI2u+RqUkDsIfbVtDyjVAiMCBoiJ8BFeC4aQTzCoxS2BBcOXNKG9Tq\/dYg9xD5\/6\/GYEwLjEcyCpsBFQQAtsjeCeFhUU3QSFw5c0kqtNMQeYi+ibUJgRNCAyL4kGhpBPEqCKlgxLhy4pA3mcEHDEHuIfcEm5Cc7xN4Pjj6sgAsfKKbbaATx0G8fNe8oCY9g8RJcOHBJW9yzMBYg9hD7MC3L0SoExhGwgMnBRUBwPzDdCOLBAv\/AAw\/QwoUL6bXXXqOnnnoq8XbR8Ii6l+DCgUtad0\/KyQGxh9iX09IySoHAiKABw\/gl0RAnHi+9\/k5Jpecr5vLm88Zl5OvHv\/GNb9CkSZNow4YN1NzcnM94BblcBNwlbQVVsSoSYg+xt2oooRNB7EMjbG8fXNhjlTdlnHis\/auf0tq\/OpjXZPB8K35tFq34tSvHlaMP5wd3wmMBLgLuktaji15NQewh9l4bVF5jEJi8yPnPBy78Y2paTIrsJUf3HNmb0T3X48knn4yqx8P57e3t4cHzVIKLgLuk9eSedzMQe4i990aVxyAEJg9qYfKAizC46lYbQTxef\/11uv\/+++m+++6L5uy3bNlCd999NzU1NYUH0EMJLhy4pPXgWhATEHuIfZCG5WoUAuOKWLj04CIctspyI4jHV77yFVqwYAF1dHRE1dq5cycdPHiwbhbpuXDgkjZ868lXAsQeYp+v5XjOBYHxDGgBc+CiAHiWWRtBPCyrKjaZCwcuaaVW2JvY85DOokWLaHh42Kmuc+bMoV27djnl8Zm4EUj0iUdVtiAwVSE\/vlxwEZ4L9DvhMc4qwYUDl7RZ5Vb1u3exX7FiRW1YJ6tSDCCv5ITYZyHV+L9DYORwDC7Cc9EI4hEepbAluHDgkjas1\/mtQ+wthvHfPPwuvXH4Xbqic0p+pJEzFQEIjJwGAi7Cc9EI4hEepbAluHDgkjas1\/mtexN73QU1pB8X5atofnBwUMQBDGkkssjvvuswvXn4PZq7bAa9+dK70f\/x8Y8ABMY\/pnktgou8yNnnawTxsK+tzJQuHLiklVlbIoh9SmT\/gwePRryxwB\/a9xY9edcI3bKpFRF+gNYMgQkAak6T4CIncA7ZGkE8HKorMqkLBy5pRVaWPIs9b71YuXJlZl37+vrEbM9II\/Fbt74YCbuK5iH4mdTmTgCByQ2d94zgwjuk4ww2gnjYomQ7mlv2SXwuHLiktcWl7HSlR\/ZlVzCrvDQSV894nm7fddWYSJ6jfR7e\/8XuaYjws8B1+B0C4wBW4KTgIjDARNQI4mGLEsTeFqmw6YKIfViX\/VpPeuie33mcWNjvfOaT4wrk357f8Tp9ftfVfp2ZwNYgMHLIBxfhuWgEsU+74pYvyOnt7aXR0VHikVyur1qnperOKM+bN4\/WrFlDTzzxRG1UeGBggLq6uoKT4MKBS9rgjucsIKjYxw3rl0WkLR6KxO\/\/yVZqbW2lyZfOirI+edfhaAV+nKBzZP+tW38SDe9f0zXNtiikS0EAAiOneYCL8FzEicfpY3IvwWFEVN+o0Em64lYt0O7u7o5Em18KlNjzsbr8ErB+\/XriM1b6+\/uppaUlmtbFMH7YdhdM7Fnod+zYUXub42qYjSBs1eysv\/TC07RkcQ9t+PgZ4oftY1\/7f9T0X24knq\/nT1L0zoL\/zev+hW7Z1AbBt4M6NRUExgOInkyAC09AppiJE\/vjf\/w1Ov6d+8MXnrOEabfdR9N+82tjcsddccvfrV69unblrf73Cy+8QNu3b4+ieT5DXx\/if\/TRR2nWrFmlRPVcCZdo3SVtTniDZwsi9vW29W72KwvoZ\/tepY+0\/Xc68deP0+W\/\/1Pa\/OvvRCKettVORfhZ6YKz2AAFQGDkkAguwnORFNm\/94rc6P7DH501LrpnpMyI3Jyj18V+z549tG\/fvprY679B7MO2O4j900\/TlH9ZQJe8NoVmLvlHeul\/XEnvX91D29Z1jVucF0cFr9Dnuf2pbedi0V6BtgqBKQCe56zgwjOgMeYaIVJU0bF5xS0i+\/DtJ08JQcSeHalyGF9fK5C1RoAfOhb7los+Rc2f20U8lDbyR4\/S7r\/7Syux57qah+9gHt+9KUJg3DELlQNchEL2P+w2gtgnXXHLteS5+M7OTszZh29K1iUEE3sl+Oa++yzxtfY8ISG\/VS5fvpzWrVsXpVD\/bm9vj83BD92b\/+c2+uy8brrglx6M5u2f6bmNvv\/8YLQS\/+K2c61dUiv4L277cLQt75ruaU75rQtqwIQQGDmkgovwXDSC2Kddcauvxv+d3\/kd4r\/53vvm5ubaXDmjrFbj8\/y9CtJCa4Ri14UDl7ThW0++EoKKfT6XiuXiBqPPCWWt8GQSz3vqVpp0ZSv9+7V7osL3ffl36bUXb6W\/+sJJWvFrV9Lf\/Ovx6PvPfnwatU0774N\/T411lKP853ccJxZ+\/jdH+Zd\/ZgqpF4BitWvc3BAYOdyCi\/BcNIJ4hEcpbAkuHLikDet1fusNJ\/Ys7vzhrRz8Mf82oWISW\/\/sJjq\/8xL61ac3Udu0Jup4sY3O\/8ufUtM9r9Lsz95Mh4+fpL\/51zdo\/6G3ouynT5+mKy89eylOx6wp9Jmrp9H+Qyeo5aLJdOWlF9aKeGnfW3TR2x+K5vQvfntS9D2LvxotaLn+I\/TmS+\/R9Kuaanl47l\/\/XPCxs\/ka\/QOBkcMwuAjPxTPPPEM9PT20d+\/eaMsvPuUjUNt2\/f3vx3Jw5MiRmlMjIyN1z5c3sXe5z573V4a6CMeM5DnSP3jwYOLxvP\/8tz+hk7+\/kK66bYR+\/fc\/Scde+DHdcOou+swFF9EvL\/oLOuc3HxzTCkf\/\/XT098snTtMzI+9E\/\/\/0ZefR\/\/67N6L\/P\/tv70S\/z7\/6Lfr7H79F\/7X1Rfr7Vy6m\/\/zz43TdC0dp8ozzaNKhj9H7b52mN6f8nN5+p4VOvPbR6LsLrh2O\/n\/hJa\/Q6I9\/gS6c\/gq9f+5MmvTukSgNf\/g7\/vDf\/G\/1PX\/n+mJw4fSjdOI1ORf7\/Pz0afrQ5MnlP\/UocRwC4CJso3j5Q\/9Ef3x8FcQ+LMyp1pXYL7r4MZpy5mz\/qn90sefvd57\/W3XNlzexN4FKW6Dncue9a1vII\/bf+423aP7Su6nps79Kx865np76n+\/SjI+cpl++ZhFduOiRaN991ufUS9+l9157mt7\/2b8RvXaaJp3+Eb3\/9tkXA3p\/Jk268FV6\/8QldN4vfDJaF8DfHX5tKl31n6bSpAsmR2lZ+PklgAVcHWChvlfl62Kv+3Ti1fGNNctnfqnIky\/Lbp7febTkxFsn6MIpF9JkCH4eCL3lARfeoEw09C+jz9HA9+6sa\/EIj1LYEpTYPzD\/O9R6Wfroysi\/jdDdT9xW13wFEfsq99m7DuNzZK\/E\/ur5v0JNn\/hf0YE6PJz+K30v0lt\/vYWm3LhwnOC\/\/7MROnX4u3Tq0A\/ptNobO+lIJOIfOu+6qJWee8VnaNL5rZS0PzVsU64v6xg6lsMXuAjPRSPMAYdHKWwJLhy4pA3rdX7rDSf25rB91gI9Jfa\/ce8j0Qp6XpHPJ+Opg3KOfbO3hu7kS6+gc6\/8OJ06\/Cc06aJXI5HnaJ3F\/fzrbrcaAchPVWPnhMDI4RdchOeiEcSDg7olS5bQqlWrojlv3m63YMEC6ujosAbQ3JNvndFDQhcOXNJ6cC2IiSBiz55Wtc\/edeudEvvbNv4DTT\/\/a9Q8\/6fEt93NXTqjdnoeD7u\/+RcP0M\/feSYSe5p0hJo+8ZVI5M3zooOwNAGMQmDkkAwuwnPRCOIBsQ\/fTnyWEEzs2UnVoHWHt23b5vTml6eyLofqKLHnq2wvfPUXo8h+w3Wzx515\/\/Y\/LItcmXT+ZdFQPz5+EYDA+MWziDVwUQQ9u7yNIPa8z3737t3RRTaPPPIIPfbYY\/SpT30qWnzNt93p++X1ffdqgTYjtWjRIhoeHo4uxeF8fHY+71JQn5B64cKBS1q7FlB+qqBiX3513EvUxX7amS9Fi+O2fGlR7fS8kz\/+PXrv1f2R4Sm\/9GA0B4+PfwQgMP4xzWsRXORFzj5fnHjwOiDJH7Pvi4vs2X++5IYFfOnSpTQ0NETTp0+PRF0tzOapVX4Z4HS8pU1dmsM34umHoJlnpvjGxkXAXdL69tOXPYj9Bwv0OLJvmf1PdOT\/Pkh\/9JX7iIf1Wz7xT\/Tzn41EQ\/YfvsR+HsoXORPJDgRGDtvgIjwXceLBgcXJH28MX3jOErgf1Ec148ReHZGr\/8YizgKvtlvr8\/T8m35Dnu4aY6TfkJfT7cRsLgLukta3n77sBRN7fdjGdDbkPntXYPTInhfo8cU23771J7Tw0UH6aOcC+tD5rYjmXUHNkR4CkwO0QFnARSBgNbNJkT0HF1I\/Zl+YNmdvir0+NM\/146F\/jvr5o4s9vxRs3ry5BoF+nK5vXFwE3CWtbz992Qsi9idPnqxdhDB\/\/vzaKs3Zs2ePGc7xVYkidpLE3vVc\/CI+IC8RBEZOKwAX4bloBPFwEfukCN28614fAUBk77cdBhF7c5+9vv0tNIGu8Jhiz2faP3nXYVp19BpXU0hfAAEITAHwPGcFF54BjTE3kcTenLPXd2rpw\/i8OE+JPV+Mw1v5+MNz+\/y3748LBy5pffvpy14pYq\/vfWfQ9Lc3XxXJawdinxc5v\/kgMH7xLGINXBRBzy5vI4iHGsF99tlna6vx1T57PernG0f1aV01hM\/fq8CQUXv44Yej20rVCn\/ev\/+d73yHNmzYEN2W5\/vjwoFLWt9++rIXROzZOf0kO13g9+zZM+ZWOl8VyWvHFHu+tIajex7Gx6c8BCAw5WGdVRK4yEKo+O+NIB7FUajWggsHLmmrrVVy6cHEXp+37+rqisSfF17ob3USQIHYS2ABc\/YyWDjrBcQ+PBuNIB7hUQpbggsHLmnDep3fejCxz+9SuTlNsef5+jcOv0uf33V1uY5M8NIgMHIaALgIz0UjiEd4lMKW4MKBS9qwXue3HkTsq7wIxxUKiL0rYmHSQ2DC4JrHKrjIg5pbnkYQD7cay0vtwoFLWnk1PesRxF47VIf32fONd\/xBZF9uk4XAlIt3WmngIjwXjSAe4VEKW4ILBy5pw3qd37pXsdfPpE9zqa+vLzo6UcLHjOzV9ba3bGqT4N6E8QECI4dqcBGei0YQjyovwrE5StfcEWCy6sKBS9rwrSdfCV7FXrmQNoyfz81wuSD24bB1sQyBcUErbFpwERZftt4I4gGxD99OfJYQROx9Ohjalin2+l32ocuG\/f9AAAIjpzWAi\/BcNILYl33rncKMd3TdeOONdOLEiejAHf7wATy8P58\/fFMeH8muvkvaAebCgUva8K0nXwkQe2POHmKfryEVzQWBKYqgv\/zgwh+WSZbixOPNw++GL7hACRe3nTsmd5m33vHteL29vbR+\/fqakLMzLPabNm2q3aJn3ra3ZMkS4sN5+AAf8+Mi4C5pC0AcNGtQsdfn8Plt69ChQ6IO1GFkEdkHbV\/WxiEw1lAFTwgugkMcO4zPB3r9YP3R8IXnLGHu0hk0d9mMWu4yb71jEdfP12fx5b\/vvfde+vrXv07q5D79fJebbrqJIPb\/QXYwsVd3FvP9xIsXL44W5KmhFR5WkbpAb\/WM54kX513TNS3nI4FseRCAwORBLUwecBEGV91qUmTPZ3xI\/UxtO5f06N7lIpyit9498cQTYwJFJfa6vnR0nL2GXN3FArEf25KCiL2+QM+86U4\/OjfEeceuD4oZ2UPsXRH0kx4C4wdHH1bAhQ8U0200wrCwi9gXvfUOkX3xNgmx1+bs+c2V5+wR2RdvWK4WIDCuiIVLDy7CYassTySx93HrHQ\/PL1q0aMwIMWOJOXv7thpE7Ll4tQ9SH2ZRUX53dzfxefkSPnpkr8T+9l1XER+wg095CEBgysM6qyRwkYVQ8d8bQezLvvVOYcbo33PPPfTiiy\/S3XffHZFhrsbnIX3dv6GhoXGL9Fw4cElbvHWEsRBM7NldnRzl\/sDAgBihZ58g9mEalqtVCIwrYuHSg4tw2DZSZB8epbAluAi4S9rfMla\/AAAgAElEQVSwXue3HlTs87tVXk6IfXlYp5UEgZHBA3sBLsJz0QjiER6lsCW4cOCSNqzX+a1D7LU5e4bx27f+JLrL3txTmh9i5LRBAAJjg1I5acBFeJwbQTzCoxS2BBcOXNKG9Tq\/9WBif+DAgegQhNHR0XHe8Ra8wcFBkrYaH2KfvyEVzQmBKYqgv\/zgwh+WSZaUePD5I62treELRAnjEOCDenhL4N69ezM5gNgnNCC1MELSfvqktq4P40Psq+sRIDDVYW+WDC7Cc8FCs2zZMtq\/f3\/4wlBCIgLXX399dDhP1gdin4BQvV6Eo8R+1dFrsrjH754RgMB4BrSAOXBRADyHrCz4\/F\/Sh3l444036JJLLqHJkyc7WEZSWwR4VMVmZAVinxHZqyMMbYGvIp0e2b95+D168q7DBLEvnwkITPmYp4nMyy+\/TG1tbRCZCmnBM1Eh+EbREPsULqSdlJfkKsRexgOFjk0GD+wFuJDBBXiQwQN7AbHXuFBD93ysYdZH6gI9juz5MgpejY9PuQigYysX77TSwIUMLsCDDB4g9nJ4KOSJHtm\/9MO36fmdxyH2hRDNlxkdWz7cQuQCFyFQdbcJHtwxC5UDkX0CsmkL9PIO7+vX5aadwqd2AuzevbvmXV9fX+ItexD7UI+Hm110bG54hUwNLkKia28bPNhjFTolxL4ksec9+3zG\/rp166IS1b\/b29vHeaDfxBT3u5kBYh\/6MbGzj47NDqcyUoGLMlDOLgM8ZGNUVgqIvYG0Hn2nkZAWacflU5fq8A1HTU1NtfuK4y7T4ReD1atX04YNG6wO7YHYl\/W4pJeDjk0GD+wFuJDBBXiQwQN7AbHPEdnnoW\/t2rVRthUrVkT\/N\/\/WbbpOE+hi\/487jtMbh9+lz++6Oo+byFMAAXRsBcDznBVceAY0pznwkBO4ANkg9gFAjTPJ4j5r1qzabXkc6R88eDB2Ht4cXcha+a\/E\/saHPkKv\/8359M9\/+1Na8redJdUMxSgE0LHJaQvgQgYX4KFaHo4cOVJzwOVo3Wq9Ti492Nn4PivsIvacls\/j14f89b9Nv5TY\/\/l5X6XZp3+Zppz5KE3\/rRdo4cKFPqsAWxkInDp1io4dO0YzZ87EQS4VtxZwUTEBHxQPHqrlYcuWLbR169YxTtico1+t13Uk9npkrqLyRx99NKqBzTC+WVV9cV\/cgj09suc99pe1Xkaf\/d0PRaKDT3kI8C6Ko0ePRkdX4mjQ8nCPKwlcVIu\/Kh08VMsDR\/Yquuc7DDZu3Gh1aU61XteR2Me5ag7bm5F+GrhZC\/b0OXsW+6lt59Itm9qk8tWwfmHIUg614EIGF+BBBg\/sBebsS+LCduud2mPf2dkZze\/b3L4HsS+JxIxi0LHJ4IG9ABcyuAAPMniA2OfkQb0hcXa+AndoaIhs9sMnHaqjBF1dumMeqjNv3rza\/H2cyxD7nER6zoaOzTOgBcyBiwLgecwKHjyCWdAUIvuCAErIDrGXwAKiSRksnPUCIiODDfAggwdE9nJ4KOSJLvZP3jVC13RNo7nLZhSyiczuCKBjc8csVA5wEQpZN7vgwQ2vkKkR2YdEtyTbEPuSgM4oBh2bDB4Q2YMHOQjI8QRin8IFL6rr7e2N9rybn6yDbsqkGGJfJtrJZUHsZfAAsQcPchCQ4wnEPoELfVX8\/Pnzqb+\/n3gB3ezZs2nRokXRfvmOjg4RTELsRdCAeWIZNERe4MVLBhngQQYP7AXEPoEL84pbfV88g7Z9+\/bUFfJlUqyL\/bdv\/QnNXToDc\/ZlEvBBWejYKgA9oUhwIYML8CCDB4h9Cg+m2OuH4rheVBOablPs+UAdXqSHT7kIoGMrF++00sCFDC7AgwweIPYZPOg30+kCv2fPHtq3b5\/YyB5iX80Dho6tGtzjSgUXMrgADzJ4gNhn8GCeZsfiv3nzZqeDdMqgGpF9GShnl4GOLRujslKAi7KQTi8HPMjgAWIvh4dCnkDsC8HnLTM6Nm9QFjYELgpD6MUAePACoxcjWKCXAKM5Z68nkzpnz8P3T951mG7fdRVd0TnFSwOBEXsE0LHZYxU6JbgIjbCdffBgh1MZqSD2EPsy2tmEKAMdmxyawYUMLsCDDB4wjB\/Dg35ZTRpNfX19tbvpq6ZTDeMjsq+WCXRs1eKvlw4uZHABHmTwALFP4SFtGF8OfWc9gdjLYAQdmwwe2AtwIYML8CCDB4i9HB4KeQKxLwSft8zo2LxBWdgQuCgMoRcD4MELjF6MYM4+A0Z9WH\/btm106NAhUXvsEdl7eQ68GEHH5gVGL0bAhRcYCxsBD4Uh9GYAYp8CJe+r50twli9fTosXL47m6PkCHD4nv6WlRdycPR+T+4P1R+nOZz5JF7ed662RwJAdAujY7HAqIxW4KAPl7DLAQzZGZaWA2Ccgrc\/Zm5ffSN16B7Ev67GJLwcdW7X466WDCxlcgAcZPLAXEHuIvZzWWOeeoGOTQyC4kMEFeJDBA8Q+gweer+cz8PVhfBXld3d3U1dXlwgm1QI9RPbV0oGOrVr8EdnLwV95gmdCDieI7DO4UADpyQYGBsQIPfsFsZfxQKFjk8EDewEuZHABHmTwgMheDg+FPFFiz9faPr\/zOK06ek0he8icDwF0bPlwC5ELXIRA1d0meHDHLFQORPahkC3RLsS+RLBTikLHJoMHRPbgQQ4CcjyB2KdwceDAAert7Y2235kf3oI3ODhIzc3NlbMJsa+cgsgBiL0MHsAFeJCDgBxPIPYJXKi77CXtp09qNhB7GQ8UxF4GDxB78CAHATmeQOwTuKjHs\/ExZ1\/tgwWxrxZ\/vXRwIYML8CCDB\/YCYp8R2S9YsIA6OjrkMBbjiR7ZH9r3dnSCHj7lI4COrXzMk0oEFzK4AA8yeIDYZ\/Ag7aS8rGH8yzsvoDcPvwexr+j5QsdWEfAxxYILGVyABxk8QOwNHtTQ\/fDwcCZDEhfoQewzaQuaAB1bUHidjIMLJ7iCJQYPwaB1NoxhfGfI5GVQw\/gQ+2q5QcdWLf566eBCBhfgQQYPiOzl8FDIE4h9Ifi8ZUbH5g3KwobARWEIvRgAD15g9GIEkX0CjGmr8cucy+drdmfNmpV6PC\/E3suzUNgIOrbCEHozAC68QVnIEHgoBJ\/XzBB7wWLPQr9582bKOosfYu\/1mchtDB1bbui8ZwQX3iHNZRA85IItSCaIvQEr33S3cuXKTLD7+vpoxYoVmenyJFCjClOnTo2y33zzzVaR\/cVt59LFbR+mz++6Ok+xyFMQAXRsBQH0mB1ceASzgCnwUAA8z1kh9jkie88cjDPHYn\/8+HHi0\/v6+\/ups7PTWuzfOucVun3XVTRz5szQbsK+gQA6NjlNAlzI4AI8VMvDkSNHag6MjIxQT08P7d27l1pbW6t1LGfp55w5c+ZMzryis6kje23F\/sQ5r9Bbk16hPz\/vq3THHXfQwoULRdev0Zw7deoUHTt2LHrRmjx5cqNVr67qAy5k0AUequVhy5YttHXr1jFOQOyr5SS2dFexV8P4H19xVnAQ3ZdLKvN19OjR6K0ZYl8u9mZp4KJa\/FXp4KFaHjiyV9H9\/v37aePGjYjsq6LEvFlv27ZtteN584o95uyrYRNDltXgHlcquJDBBXiQwQN7gTl7OVyM8wRiL5icGNfQscnhC1zI4AI8yOABYi+HBy\/D+GyET9FDZF8NsejYqsEdkb0c3E1P8EzI4QaRfQ4uFGiclVfMDw0NUXt7ew5L6VlcI3u2xtfc3rKpzbsvMJiNADq2bIzKSgEuykI6vRzwIIMHRPZyeCjkiTpUB2JfCMbCmdGxFYbQmwFw4Q3KQobAQyH4vGZGZO8VzmqMQeyrwR1DljJwj\/MCIiODG\/AggwdE9ik8SDkb36apQOxtUAqfBh1beIxtSwAXtkiFTQcewuLrYh2RfQJaEHuXZoS0jAA6NjntAFzI4AI8yOABkX0MDxLOxndtHojsXRELkx4dWxhc81gFF3lQ858HPPjHNK9FRPY5Ivu8YIfKp4v93KUzaO6yGaGKgt0UBNCxyWke4EIGF+BBBg+I7OXwUMgTiH0h+LxlRsfmDcrChsBFYQi9GAAPXmD0YgSRfQaMccP6WffLe2HGwQjE3gGsgEnRsQUE19E0uHAELFBy8BAI2BxmIfYpoLHQ79ixgwYHB6m5uTlKqRbudXd3p147m4OL3Fkg9rmh85oRHZtXOAsZAxeF4POWGTx4g7KwIYh9AoT1uhofc\/aFn4ncBtCx5YbOe0Zw4R3SXAbBQy7YgmSC2EPsgzSsiWgUHZsc1sGFDC7Agwwe2AuIPYbx5bTGOvcEHZscAsGFDC7AgwweIPYWPGCBngVISBIhgI5NTkMAFzK4AA8yeIDYy+GhkCf6Aj2+8Y5vvsOnfATQsZWPeVKJ4EIGF+BBBg8Qezk8FPIEYl8IPm+Z0bF5g7KwIXBRGEIvBsCDFxi9GMGcvQajWoE\/PDycCe6cOXPGbMnLzBAwAcQ+ILgOptGxOYAVOCm4CAywpXnwYAlUCckg9ikgp+2zX7FiBXV0dJRAUXYREPtsjMpIgY6tDJTtygAXdjiFTgUeQiNsbx9in4BVve6zx5y9feP3nRIdm29E89sDF\/mx85kTPPhEs5gtiD3EvlgLQu4aAujY5DQGcCGDC\/Aggwf2AmKfcxhf6nG5iOyre7jQsVWHvVkyuJDBBXiQwQPE3oIH9TakJ922bZuY+Xr2C3P2FkSWkAQdWwkgWxYBLiyBCpwMPAQG2ME8IvsEsE6ePEn8n7oAxwHT0pPqYn\/7rqvois4ppfuAAnGojqQ2AJGRwQZ4kMEDIvsUHrJut9uyZQvNmzdPxMsAxF7GA4WOTQYP7AW4kMEFeJDBA8Q+gwc17MGivmbNGmpqaqIDBw5Qb28vXXrppSL32SOyr+7hQsdWHfZmyeBCBhfgQQYPEHsLHlSEf+zYsUjg+cCdvr4+4n32Uj6I7GUwgY5NBg+I7MGDHATkeII5ewsuVDQ\/OjpKAwMD1NXVZZGrvCQQ+\/KwTisJYi+DB4g9eJCDgBxPIPYZXKxdu5Y2b94cRfM33HAD9fT0RHP1alhfApUQewksYJ5YBgtnvcCLlww2wIMMHjCMn8KDPnw\/NDRE7e3tUWr1Pf97cHAQC\/TktOXKPUHHVjkFNQfAhQwuwIMMHiD2GWK\/e\/duWrhwYWwqrMaX04ileIKOTQoTiOylMIFnQgoTOEFPDhMFPMEwfgHwPGZFx+YRzIKmwEVBAD1lBw+egPRgBnP2HkB0NcHrAGbNmpW40I8P8+nv7yceWVCftB0Autjf+cwn6eK2c11dQnoPCKBj8wCiJxPgwhOQBc2Ah4IAeswOsfcIpo0pteAvbVU\/rwtYsmQJrVq1qrZWIM02xN4G+fBp0LGFx9i2BHBhi1TYdOAhLL4u1iH2LmgVSKsW9k2dOjWycvPNNydG9rzVb\/Xq1bRhwwarBYAQ+wLEeMyKjs0jmAVNgYuCAHrKDh48AenBDMTeA4g2Jljsjx8\/Ti0tLdEQfWdnZ6LYMyk8AmC72h9ib8NA+DTo2MJjbFsCuLBFKmw68BAWXxfrEHsXtDykVfPxaWK\/c+dOWrlyZa20OXPmpAq\/Lvb\/7c+mRHP2M2fO9OAtTLgggI7NBa2wacFFWHxtrYMHW6TCpDty5EjN8MjISHROzN69e6m1tTVMgYGtnnPmzJkzgcsYY16\/9pYjdX0ffpYfNmLPUT2f1qcO7jH\/NsvQxX7H+b9Fb53zCt1xxx2J2wazfMTv+RA4deoU8bHK\/KI1efLkfEaQywsC4MILjIWNgIfCEBYywFvEt27dOsYGxL4QpGMz65G5GZXbiL3pCs\/hL1++nNatWxe7YE8X+08\/\/CZd3PbhSHAQ3Xsk1cIUc3v06NHorRlibwFYwCTgIiC4DqbBgwNYAZJyZK+i+\/3799PGjRsR2Zs4h7rPPq\/Ypy3Yw5x9gKckh0kMWeYALVAWcBEIWEez4MERsIDJMWefAG6o++yzxN78Xf3N0wVJN+3pYr\/q6DUBmwtMpyGAjk1O+wAXMrgADzJ4YC8g9ilchLjPPk7s1XcLFiygjo4OMg\/Vybp4B2Iv44FCxyaDB\/YCXMjgAjzI4AFib8FDvd1nj8jegtRASdCxBQI2h1lwkQO0AFnAQwBQc5pEZG8BXD3dZw+xtyA0UBJ0bIGAzWEWXOQALUAW8BAA1JwmIfYZwNXbffYQ+5xPgods6Ng8gOjJBLjwBGRBM+ChIIAes0PsE8Cs1\/vsIfYenw5HU+jYHAELmBxcBATXwTR4cAArcFKIfYrY1+N99hD7wE9Minl0bNVhb5YMLmRwAR5k8MBeQOw1LlQ0z1vceFW8zUedY79r1y6b5EHSYDV+EFidjaJjc4YsWAZwEQxaJ8PgwQmuoIkh9hD7oA1sIhlHxyaHbXAhgwvwIIMHRPYGDyqyHx4edmKIj8SVENnzBTh3PvNJJ9+R2B8C6Nj8YVnUErgoiqCf\/ODBD44+rCCy94FixTbUMD7Evloi0LFVi79eOriQwQV4kMEDIns5PBTyBGJfCD5vmdGxeYOysCFwURhCLwbAgxcYvRhBZO8FxmqNQOyrxV+Vjo5NBg\/sBbiQwQV4kMEDIns5PBTyBGJfCD5vmdGxeYOysCFwURhCLwbAgxcYvRhBZO8FxmqNQOyrxR+RvQz8dS8gMjI4AQ8yeEBkn8GDvu9+9uzZtGjRIuKV+rz6fnBwkJqbm0UwCbEXQQOGjmXQEHkBkZFBBniQwQPEPoMHPhefP3zIzs6dO2nHjh2RyO\/Zs4cOHjyYeL982fRC7MtGPL48dGwyeIDYgwc5CMjxBMP4CVzoUT1H8v39\/dTS0hIJvDo1T0p0D7GX8UBB7GXwALEHD3IQkOMJxN5C7NUQfnd3N3V1dUHs5bRfUZ5A7OXQAS5kcAEeZPCAYfwUHk6ePBlF852dnXTFFVfQ0qVLaWhoiNrb20kf3pdAJSJ7CSxgnlgGC2e9gMjIYAM8yOABYp\/Bw4EDB6i3t5dGR0epr68vGsJnoee\/16xZQ01NTSKYVGJ\/eecF9PldV4vwaSI6gY5NDuvgQgYX4EEGDxB7OTwU8gRiXwg+b5nRsXmDsrAhcFEYQi8GwIMXGL0YwZy9FxirNQKxrxZ\/VTo6Nhk8YBgfPMhBQI4nEPsMLnjL3cqVK6NU27Zto0OHDtG+ffswjC+nDYvxBGIvhgrM2QuhAs+EECKIooXlPT09tHfvXmptbZXjmIMn55w5c+aMQ3rrpGp+fvny5bR48eJozt7chmdtLGBCRPYBwXUwjY7NAazAScFFYIAtzYMHS6BKSAaxTwA57vQ8FvuOjg6xW++wQK+EJyalCHRs1eKvlw4uZHABHmTwwF5A7CH2clpjnXuCjk0OgeBCBhfgQQYPEPsMHni+nufn9WF884AdCVRiGF8CC9jbLYOFs15AZGSwAR5k8ACxt+BBDX3oSQcGBqKT9KR8IPYymEDHJoMHiD14kIOAHE8wjC+Hi9yeQOxzQ+c1I8TeK5yFjIGLQvB5ywwevEFZ2BDEvjCE1RuA2FfPAaJJGRwoLyAyMvgADzJ4wDB+Cg\/6anxega9\/pN56d03XNLplU5uc1jXBPEHHJodwcCGDC\/AggweIPcReTktsAE\/QsckhEVzI4AI8yOABYh\/Dg35iXhpN6mKcEFSqUYXh4eHI\/Lx581JP7FPD+IjsQ7BhbxMdmz1WoVOCi9AI29kHD3Y4lZEKc\/YJKKcN44ckRr9al1f8q79bWlqiE\/ziPhD7kIzY20bHZo9V6JTgIjTCdvbBgx1OZaSC2JeBcsEy1H7\/pGt1IfYFAfaUHR2bJyA9mAEXHkD0YAI8eADRkwmIfQqQ+n32ZjI+I39wcJCam5s9UZFsBmIfHGIvBaBj8wKjFyPgwguMhY2Ah8IQejMAsU+AUh9Onz9\/PvX399OCBQtInaCnzsn3xkSCITWd0N3dnXiQj4rsmz\/7Nt2y8exq\/JkzZ4Z2DfYNBNCxyWkS4EIGF+ChWh6OHDlSc2BkZAS33sXRYc7Z8w14s2bNigSX35C2b98e\/Jpb9cLB\/iUN4fNvSuwPTP4+PfWRTVF17rjjDlq4cGG1LW2ClX7q1Ck6duxY9KI1efLkCVZ7WdUFFzL4AA\/V8rBlyxbaunXrGCdwxa3BiSn2PJR+8ODBaJGcz3325lTBtm3bopv1bIVeF3uO7D+17JxaZI\/ovtwHjTk7evRodFc0xL5c7M3SwEW1+KvSwUO1PHBkr6L7\/fv308aNG3GffRwlHM3zxxT4PXv2RBfkpEXbRSi2WYGv28cCvSJo+8uLIUt\/WBa1BC6KIugnP3jwg6MPK5izT0HR3AbH4r9582bibXBDQ0PU3t7ug4NxNric0dFR65cJJfZzl86guctmBPEJRrMRQMeWjVFZKcBFWUinlwMeZPDAXkDs5XAReWIeqKPcS1v9D7GXQSI6Nhk8sBfgQgYX4EEGDxD7FB7q8Wx8RPbVPljo2KrFXy8dXMjgAjzI4AFiD7GX0xIbwBN0bHJIBBcyuAAPMniA2MfwIOFsfNfmgWF8V8TCpEfHFgbXPFbBRR7U\/OcBD\/4xzWsRc\/YJyFV1Nn4eIiH2eVDznwcdm39M81oEF3mR85sPPPjFs4g1iH0R9ITkhdjLIAIdmwwe2AtwIYML8CCDBwzjJ\/DAQ\/kPPfRQbXuduUJ+YGAg8ejaKqiF2FeB+vgy0bHJ4AFiDx7kICDHE0T2BhcMyNKlS2tCbx5wY3NWfdn0QuzLRjy+PIi9DB4g9uBBDgJyPIHYa1yYh+jwT+o42\/Xr10fH2PIn6xa6sumF2JeNOMReBuLJXuDFSwZD4EEGDxjGN3iIW5RnRvoKND7lrqwrbrOaC8Q+C6FyfkfHVg7ONqWACxuUwqcBD+Exti0Bkb2GVJzYx0XxPi\/CsSUqLR3E3geKxW2gYyuOoS8L4MIXksXsgIdi+PnMDbHX0FTD+HxvvX7zHJ+Fz5fhqI\/r2fU+CYuzpcT+lk1tdE3XtNDFwX4CAujY5DQNcCGDC\/Aggwc1It3T04Nb7xQleiQ\/MjJCvb29pM\/Xq7cjdRWtBCoh9hJYwHYvGSyc9QIiI4MN8CCDB4h9Ag\/qdjv+WW2zU1H\/s88+G\/TGuzxNA2KfBzX\/edCx+cc0r0VwkRc5v\/nAg188i1jDMH4R9ITkhdjLIAIdmwweENmDBzkIyPEEYi+Hi9yeQOxzQ+c1I8TeK5yFjIGLQvB5ywwevEFZ2BDEvjCE1RuA2FfPAaJJGRwoLyAyMvgADzJ4YC8g9nK4yO0JxD43dF4zomPzCmchY+CiEHzeMoMHb1AWNgSxLwxh9QYg9tVzgMheBgeI7MGDLATkeAOxl8NFbk8g9rmh85oRUYxXOAsZAxeF4POWGTx4g7KwIYh9YQirNwCxr54DRPYyOEBkDx5kISDHG4i9HC5yewKxzw2d14yIYrzCWcgYuCgEn7fM4MEblIUNQewLQ1i9AYh99RwgspfBASJ78CALATneQOzlcJHbEyX2t++6iq7onJLbDjIWQwBRTDH8fOYGFz7RzG8LPOTHzndOiL1vRCuwB7GvAPSYItGxyeABoyzgQQ4CcjyB2MvhIrcnEPvc0HnNCLH3CmchY+CiEHzeMoMHb1AWNgSxLwxh9QaYxG9\/7ie09LvzMYxfIR3o2CoE3ygaXMjgAjzI4IG9gNjL4SK3J41AYu7KC8qIjk0OGeBCBhfgQQYPEHs5PBTyBGJfCD5vmdGxeYOysCFwURhCLwbAgxcYvRhpBJ0458yZM2e8oFGnRhqBxDqFfozb6NjksAguZHABHmTwgMheDg+FPIHYF4LPW2Z0bN6gLGwIXBSG0IsB8OAFRi9GGkEnENk\/\/TT19PTQ3r17qbW11UvDgBF3BNCxuWMWKge4CIWsm13w4IZXyNQQ+5DoarZff\/11WrRoEQ0PD0ffzps3j9asWUNNTU3jPDh58iT19\/fT7t27a7\/19fXRihUrYr1tBBJLoiFoMQcPHqTHH3+cvvjFL+KlKyjS2cbBRTZGZaQAD2WgbFdGI+iE+MheiXdnZyd1dXWR+rulpSVWwPnFYMmSJbRq1Spqb2\/PZLIRSMysZB0kAA9ySAIXMrgADzJ4YC8agQvxYh9H986dO2nfvn2x0f2BAwdo9erVtGHDBmpubs5sLY1AYmYl6yABeJBDEriQwQV4kMEDxL5CHtLEnh+QtWvX0uDgoJPYb9u2DcPHFXI6MjISrZ0ADxWS8EHR4KJ6DtgD8CCDB52Lel7bVXeRvZq\/7+7ujob1zQ+\/CKxcubL29Zw5c1KFnx+oZcuW0f79++W0LHgCBIAAEAACohC4\/vrrafv27aJ8cnGmrsRezddzBZMW6HFUPzo6Wvvd\/DsOHBZ8\/g8fIAAEgAAQAAJxCPBurXresSVO7PXIXI\/KbYQ+jiCew1++fDmtW7fOasEemjkQAAJAAAgAgUZDQJzYxwGctQI\/jRTXBXuNRjDqAwSAABAAAkCgLsTeZiieqXTdpgf6gQAQAAJAAAhMBATEi715oI4iRQ3x88E6fIjOggULqKOjoyb46lCdtAN4JgLBqCMQAAJAAAgAAfFiD4qAABAAAkAACACBYghA7Ivhh9xAAAgAASAABMQjALEXTxEcBAJAAAgAASBQDIEJL\/a8+G\/z5s0Rinx6G8\/741MuAnHrMgYGBmIPTSrXs4lTGj8Hs2bNGoO5zkvW4VQTB6nwNY3jQh2dq0rnu0GGhoawndgzHVmXrtXzMzGhxV4\/WveFF15wOmbXcxub0OawPbJa+tULr\/mCxd\/zh2+M1P9drbeNXXoSF3z+CN+Cl3R7Z2OjUk7tbHZz1fMzMaHFXidOEa1W9ZfTvFAKI8AvXXwMZdKpiEApDAIqSpk6dWpUwM0331yL7NVvLC482oUXsjAcKKtpXHCauGg\/rEewzgjo97CwRvBV6\/X6TExYsU96i1NX6aKpl4dA2hFm8jsAAARASURBVMVG5Xkx8UpigTl+\/DjxkDBvX9XbvnnyJE6iDNs+0rgw+6qwnsC6joDeN\/GR6vpprPX2TEx4sdcjebw9V\/Og6+sm2AOcjVAuD3FiYkbyLEZLliyhVatWYZ44ID1xXGBNS0DAU0ybl67V+zMBsf\/gMB4MlVXzQJlHIRc5GrmaGtR\/qRB7ORwmcdHb20vr16+vTanof8vxvnE8ibuLBWJfp\/xiGF8ucfrCyebmZrmONohnSQJTz0OW9UqN7ZA9FkyGYzjp0rV6n9qasJG9GcljgV64h8fVMhbsuSJWLH3S0LE+bI8FesUwts3tIvbmVknbMpAuGYG0kUVzKqvenokJLfbYelf9Y2+u+jbnyar3sPE9SBKYet5mVK+sxXFhjnTx30uXLsU++wAkZ126Vs\/PxIQWexXd41CdAE+Ng0lzAVJfXx\/2EzvgVzRpktjX8wEiRTGpKn8SF+ahOjgAzD9DWZeu8ZRiPT8TE17s\/TcZWAQCQAAIAAEgIAsBiL0sPuANEAACQAAIAAHvCEDsvUMKg0AACAABIAAEZCEAsZfFB7wBAkAACAABIOAdAYi9d0hhEAgAASAABICALAQg9rL4gDdAIBEB3tfLJ6eNjo4mpnnwwQfpW9\/6Fq1bt67UY21tz0aot73JaI5AoFEQgNg3CpOox4RDQMopaq7n5uO61gnXVFFhAQhA7AWQABeAQB4EpIi9q3i7vhzkwQZ5gAAQGIsAxB4tAgjUKQJxYq+f3z19+vTo\/u358+fT4OBgNPzP19kODQ3Rc889RytXroxqbt4yaB7gMjAwULvn3oTKPAGRf1cHw+zevTtKrspsb2+vZXd9QahTiuA2EBCDAMReDBVwBAi4IWAr9myVxZ5PAFPXCSsBN48nZhHesWNHLX3W8cVxlxaZful3gjc1NUWVxGVHblwjNRAoigDEviiCyA8EKkLAVuy7u7trkbkpsvrxrDwC0N\/fTwu0a5+zhNkUctsris0bxCqCEMUCgQmDAMR+wlCNijYaArZiv2LFiuge9Djh1sX+pptuiob9h4eHx0E1Z86cWrSv\/5gUtff09ETJ4obw+fu44f9G4wf1AQKSEIDYS2IDvgABBwR8i\/21114bbe1bv3597eUgy504sVd59Ll7U\/Qh9lnI4ncg4BcBiL1fPGENCJSGgG+xV5G9PuyfVZk0sTdFv7OzszadgGH8LGTxOxDwiwDE3i+esAYESkPAt9h3dXURi\/dDDz005q70tDu+zUNy4ubs44QdC\/RKayYoCAhECEDs0RCAQJ0iEELsGQoWfLUtj\/82t+bpcClx1xf1xd0Lbt6\/LuWMgDqlHm4DAWcEIPbOkCEDEAACOgKue+ZxqA7aDxAoHwGIffmYo0Qg0FAIuIq368tBQ4GFygCBihCA2FcEPIoFAo2EAC7CaSQ2UZdGRABi34isok5AAAgAASAABDQE\/j\/ILPiXk3IuEgAAAABJRU5ErkJggg==","height":305,"width":507}}
%---
%[output:9791a366]
%   data: {"dataType":"textualVariable","outputData":{"name":"xb_max","value":"0.0258"}}
%---
%[output:15857bb6]
%   data: {"dataType":"textualVariable","outputData":{"name":"theta_max","value":"0.0873"}}
%---
%[output:6adab004]
%   data: {"dataType":"textualVariable","outputData":{"name":"disp_percent","value":"27.3553"}}
%---
%[output:0e6468c6]
%   data: {"dataType":"textualVariable","outputData":{"name":"rot_percent","value":"5.5556"}}
%---
%[output:19a15f7d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"exponent":"3","name":"Ke","rows":4,"type":"double","value":[["0.4079"],["3.3863"],["0.0630"],["0.9850"]]}}
%---
%[output:973f237d]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_o","rows":4,"type":"double","value":[["0.0200"],["0"],["0.1309"],["0"]]}}
%---
%[output:63b81036]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"Ke_min","rows":3,"type":"double","value":[["6.9631"],["57.7303"],["16.6893"]]}}
%---
%[output:0ea290d6]
%   data: {"dataType":"textualVariable","outputData":{"name":"obs_rank","value":"3"}}
%---
%[output:3c41c227]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_omin","rows":3,"type":"double","value":[["-0.0020"],["0"],["0"]]}}
%---
%[output:10205beb]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('load')\" style=\"font-weight:bold\">load<\/a>\nUnable to find file or directory 'noise_data.mat'."}}
%---
