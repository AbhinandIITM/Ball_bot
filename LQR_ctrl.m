clc ; clear;
%%

% --- 1. Define Parameters ---
M_w = 4.3;   % Mass of one wheel
M_r = 10.12;   % Mass of robot body
m_b = 2.71;  % Mass of the ball
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

% --- 3. Derived Parameters ---
l = h + r; % Height of ball CoM from wheel center


M1 = [(I_b/(r^2)) + m_b, m_b*l;
      m_b*l,            I_r + m_b*(l^2)];

% M2 = [[1,  1+(g*m_b)],
%       [1-(m_b*g), 1 - (M_r*d*g)-(m_b*g*l)]]
M2 = [1,             1 + (g*m_b);
      1 + (m_b*g), 1 + (M_r*d*g) + (m_b*g*l)];

% M3 = [[m_b*l],
%       [-m_b]]
M3 = [-m_b*l;
      -m_b];

% --- 5. Calculate A and B ---
% This is the numerically preferred method for A = inv(M1) * M2
A1 = M1 \ M2;

% This is the numerically preferred method for B = inv(M1) * M3
B1 = M1 \ M3;
%%

A = [0 1 0 0;... %[output:group:229a734b] %[output:60988226]
    A1(1,1) 0 A1(1,2) 0;... %[output:60988226]
    0 0 0 1;... %[output:60988226]
    A1(2,1) 0 A1(2,2) 0] %[output:group:229a734b] %[output:60988226]
B = -[0 ; B1(1); 0 ; B1(2)] %[output:06e878d2]
C = [0 0 1 0] %[output:717ff987]
D = [0] %[output:088d95c6]
%%
[num,den] = ss2tf(A,B,C,D) %[output:89545353] %[output:7e86f660]
%%
% --- 2. Create the individual Transfer Functions (TFs) ---

% H1(s) = Y1/U (from num row 1)
H1 = tf(num(1,:), den) %[output:6b1db27e]

% % H2(s) = Y2/U (from num row 2)
% H2 = tf(num(2,:), den)
% 
% % H3(s) = Y3/U (from num row 3)
% H3 = tf(num(3,:), den)
% 
% % H4(s) = Y4/U (from num row 4)
% H4 = tf(num(4,:), den)
%%
sys = ss(A, B, C, D);
%%
rank(ctrb(A,B)) %[output:7be9b0ee]
rank(obsv(A,C)) %[output:914ed936]
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
R = 1e-3;
%%
Q = [1000 0 0 0 ; 
    0 100 0 0;
    0 0 1000 0;
    0 0 0 100];
%%
K = lqr(A, B, Q, R);
%%
% Closed-loop system
Acl = A - B*K;
sys_cl = ss(Acl, [], eye(4), []);
x0 = [0.005, 0, -deg2rad(2),0]' %[output:323ef937]
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
if ~isempty(te) %[output:group:9245fc52]
    if ie(end) == 1
        fprintf('\n⚠️ Simulation stopped at t = %.3f s: |x| exceeded %.4f m.\n', te(end), x_limit);
    elseif ie(end) == 2
        fprintf('\n⚠️ Simulation stopped at t = %.3f s: |theta| exceeded 15 degrees.\n', te(end));
    end
else
    fprintf('Simulation completed without violating limits.\n'); %[output:04f8a943]
end %[output:group:9245fc52]

figure; %[output:5b30434f]
plot(t, y); %[output:5b30434f]
xlabel('Time (s)'); %[output:5b30434f]
ylabel('States [x, xdot, theta, thetadot]'); %[output:5b30434f]
legend('x', 'xdot', 'theta', 'thetadot'); %[output:5b30434f]
title('Closed-Loop Response with LQR (ODE45 + Safety Constraints)'); %[output:5b30434f]
grid on; %[output:5b30434f]

%%
xb_max = max(y(:,1)) %[output:97e24fca]
theta_max = max(y(:,3)) %[output:0eb22992]
% xb can be at max d/2 before the ball falls
disp_percent = (xb_max/(d/2))*100 %[output:35a48901]
rot_percent = (theta_max/(pi/2))*100 %[output:534cccc8]
%%
%[text] ## Observer
G = 0.07*eye(4);      % or your actual process-noise input matrix
Q = 0.02*eye(4);     % process noise covariance
R = 0.08*eye(1);     % measurement noise covariance
N = zeros(4,1);               % no correlation
[L,P,E] = lqe(A, G, C, Q, R, N);   % or simply lqe(A,G,C,Q,R)
%steady‑state Kalman gain L, the error covariance P, and a vector E
%containing the eigenvalues of A-LC
Ke = L %[output:4aab1181]
%%
x0_o = [-0.002;0;deg2rad(5);0] %[output:4ac57f12]
%%
%[text] ## Min order observer
% measured index
ia = 3;

% unmeasured indices
ib = [1,2,4];

% A partitions
A_aa = A(ia, ia);         % 1x1
A_ab = A(ia, ib);         % 1x3  (measured <- unmeasured)
A_ba = A(ib, ia);         % 3x1  (unmeasured <- measured)
A_bb = A(ib, ib);         % 3x3  (unmeasured <-> unmeasured)

% B partitions
B_a = B(ia);              % scalar (1x1)
B_b = B(ib);              % 3x1

% C partitions (measurement y = theta)
C_a = C(:, ia);           % 1x1  (should be 1)
C_b = C(:, ib);           % 1x3  (will be zeros for C = [0 0 1 0])

%%
rank(obsv(A_bb,A_ab)) %[output:724f0e11]
%%
    G_min = 0.01*eye(3);
    Q_min = 0.5 * eye(3);    % increase process noise (trust dynamics more)
    R_min = 0.005;           % decrease measurement noise (trust measurements more)
    
    [L_min_fast, P_min_fast, E_min_fast] = lqe(A_bb, G_min, A_ab, Q_min, R_min, zeros(3,1));
    Ke_min = L_min_fast;
    eig(A_bb - Ke_min*A_ab)    %[output:1fafbf59]
%%
A_hat = A_bb - Ke_min*A_ab %[output:6c24be73]
B_hat = A_hat*Ke_min + A_ba - Ke_min*A_aa %[output:720a56da]
F_hat = B_b - Ke_min*B_a %[output:407f6c67]
% measured index and unmeasured indices
ia = 3;
ib = [1,2,4];

% Ke_min is 3x1 from your reduced LQE
% Build C_hat and D_hat to map to original ordering [x; xdot; theta; thetadot]
C_hat = zeros(4,3);   % 4 rows (full state), 3 cols (x_b_hat)
C_hat(ib, :) = eye(3);    % place identity rows at positions 1,2,4

D_hat = zeros(4,1);
D_hat(ia) = 1;            % measured y in index 3
D_hat(ib) = Ke_min;       % place Ke_min into rows corresponding to ib

%%
x0_omin = [-0.002;0;0] %[output:26003d81]
%%
%[text] ## FFT
y_noise  = out.y_noise.signals.values;
ti = out.tout;
Ts = ti(2) - ti(1) %[output:59de8216]
Fs = 1/Ts; 
N = length(y);
Y = fft(y - mean(y));        % remove DC
f = (0:N-1)*(Fs/N);          % frequency axis
P = abs(Y).^2 / N;           % power spectrum
figure;  %[output:8e97ee7e]
plot(f(1:N/2), P(1:N/2)); %[output:2a3ae915] %[output:8f59b930] %[output:8e97ee7e]
xlabel('Frequency (Hz)'); %[output:8e97ee7e]
ylabel('Power'); %[output:8e97ee7e]
grid on; %[output:8e97ee7e]
title('Noise FFT'); %[output:8e97ee7e]

%%
fc = 10;
tau = 1/(2*pi*fc) %[output:861b919d]
num = [1] %[output:5ca1803f]
den = [tau   1] %[output:1ad2475a]

% 10 Hz is far above the 1.5 Hz interference
% the robot's physical angle dynamics are active around 3–6 Hz → fully preserved
% it significantly reduces high-frequency measurement noise
% it barely introduces delay (a few milliseconds)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:60988226]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["-11.8839","0","-2.5891","0"],["0","0","0","1.0000"],["47.6759","0","38.6973","0"]]}}
%---
%[output:06e878d2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["-0.8402"],["0"],["4.2977"]]}}
%---
%[output:717ff987]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":1,"type":"double","value":[["0","0","1","0"]]}}
%---
%[output:088d95c6]
%   data: {"dataType":"textualVariable","outputData":{"name":"D","value":"0"}}
%---
%[output:89545353]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":1,"type":"double","value":[["0","0","4.2977","0","11.0147"]]}}
%---
%[output:7e86f660]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"den","rows":1,"type":"double","value":[["1.0000","-0.0000","-26.8134","-0.0000","-336.4381"]]}}
%---
%[output:6b1db27e]
%   data: {"dataType":"text","outputData":{"text":"\nH1 =\n \n                    4.298 s^2 + 11.01\n  -----------------------------------------------------\n  s^4 - 7.605e-15 s^3 - 26.81 s^2 - 1.288e-13 s - 336.4\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 52 46 50 57 55 55 32 48 32 49 49 46 48 49 52 55 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 45 55 46 54 48 53 48 101 45 49 53 32 45 50 54 46 56 49 51 52 32 45 49 46 50 56 56 52 101 45 49 51 32 45 51 51 54 46 52 51 56 49 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:7be9b0ee]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:914ed936]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:323ef937]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0","rows":4,"type":"double","value":[["0.0050"],["0"],["-0.0349"],["0"]]}}
%---
%[output:04f8a943]
%   data: {"dataType":"text","outputData":{"text":"Simulation completed without violating limits.\n","truncated":false}}
%---
%[output:5b30434f]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAcoAAAETCAYAAAC2rUY7AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQuYVtV57191QEdAh0EExwEmGtDmAqnaSsZaSVtrTk8gMWnCQKx0xJSTI8c04ToYb6kCjnASVGp4lI5QBYwnJYJJm9AmikjlHLESExM0KuA43FQkoIigc5534\/pcs2Zf1t57vfv79p7\/9zw8zMy3rr+19vvf77oe19nZ2Un4gAAIgAAIgAAI+BI4DkKJngECIAACIAACwQQglOgdIAACIAACIBBCAEKJ7gECIAACIAACEEr0ARAAARAAARBIRsCpR3no0CGaPXs2rV27tlSasWPH0vz586m6utr724MPPkgtLS00b948Gj9+fLJSO4j1xhtv0OTJk72Uli5dSrW1td1SfeGFF6i5uZnOP\/\/8LnVwkL2TJBRLv8RGjRoVWC8nmRc0EZt+YRPGxBM3zm233UZLliwpJVNXV0dtbW00fPjw0t\/Stj\/Hv\/POO7ula+bNGa5YsYJGjx5dytsvDH9pllPVe8uWLb49JqifqvR1O6Gex46OjlAuabumnx0z6x+Vh15Wv7aLip\/mexs7nCb9csVVfSJuW8Qpr2LHcXTdciaUfp1YFVDvKBDKOM0WHjbMUHJM8yXFXc7FTSlI0J588kl67LHHaNasWRRX9JiWbZwoYZkyZYpXBv6kaX+VD4ufSi9O3kFCqXqGMmZRafoJJbOeOHGil5QulPrf\/WxLUK+MY3P8RMasU1TvN9OwFUqOd+utt9KkSZO6vBBF5ad\/b2uH46SZNiyXadmyZXTdddeVHKYkaSYVSv3ZtclX9TNdkJ0Ipd4x9AdZbzT19zid1qZSScPYGK68eJSmd67qtnfv3m7eQlJePTme6rOqD9v0HZOXbRxlDEwB0UVCPcBBz5JN+6u4ujGIk3eQ0VJ\/Vy9pbBuiRm50VqbI6H07qe2IE0898wMHDiyNyKg66bYt7HlQddi8ebP182fTZlHPYBw7HJWWq+9Vvy3XS7v57NrUS7XF0KFDS16lE6EMExQGtXLlylKGfp3WfOv0ewMz32BN8GYaUd\/fdddddM8993jc0g692pQ\/Kozicv3119MzzzxTGr4OG6IOMgD6A6MMoWmATEPs9ybt571Elc98ozXzUe1477330sMPP1yqp55XVFm4zUzvIoiTn9FSZeR01HCmSu\/KK68kNUzI\/WLdunXeVIHuwSxatIhuueUW708LFiyg6dOne3HCPAcboYwyluZDH6f9\/cRIN+Rx8w4SSpP3gAEDYgmlqtPFF19Mjz\/+eMmj9OsTtp5aEqHk4d2wIT7THql+zpz5xUAfalZ9O6jP+nndX\/7yl+mJJ57oMu0TJcBx7DCXM8ompX1W\/UYAmCmPzPC0AtsSfsYUa+4rPNWlD63rz7XZ56LKZ464qP7CdTfz0dvaj7MToVRAbN64zE4bNlSgCh80zGO+4ZvzIOZbrd88Sdhcno1HaVN+mzBhw2hBIhBkAEwhqK+v7zZ3zJ1F1Z3nj825ZWVYVd425fN7MDgd3aCFDdlxe3OZosoSlI8tJ70upnf27W9\/uyTeUULJHjt7Hnq\/CupPNkIZ9fZt9sc1a9b4zvf7vQjoQulXFtu8lafFL5ls7PzERDdoI0aM6CYceln8RqD4uW1oaOhSt6Ah3CDeQX1E5R0kgjbDpkHPAtfla1\/7mq9QXnLJJaXhZL3+3GcvvfTSbnFYbPfs2UP6y0yUPYpjh21sUtpnddiwYd3qrAul38unaaP9bEeULnC6HGb79u3dXnK\/853v0A033NBFjE0bxb+bouxEKOO4t6ZxN4dq2GirMEro7rjjjsCHkivll79eUQ7Dcx76Q6XipBXKOOXXvVyzjsroxSlj1ByV+Sbrlz8\/qOPGjfPEKWyoyI+X\/jf20FtbWz2R0Y2fOWxl8uK2UcJoUxY\/bzlKhHQDcu2113YRYi6r+hvXX\/cW1UhD0NArP9RKnKM8sqgy6l5y0DCVmYbp7ZrDSkEvrn4GN0oozbxthDJIBPyEUrXrjh07PC9D1U3xNYdE9Zc7vxekpELp52npL5V+i\/5MdqpPq+dJf1E1527ViJbyRPXpEtNWRnnGcexwErsV91kN6tN+eZt91+85D\/Io1fNilo8Xi5pMol42VDlM1mUVyiADbRqd119\/3ffNRK3CC3vz4YeIP+ZK2yijw0bmi1\/8Yuiq16ChEL38d999N\/EwoylCZh2ffvrpwDIGzTWGCaX+xhwWzpw71t\/y9FWWYUPmXD4\/geG0TAOnDGzY3FPQkAmv+AxbHBI0FKe3tRoqramp8ap6yimn0FVXXUXXXHONN8w1c+ZMmjp1qvddlFDq7RI1LOZCKE2WYUIZNmzoJ4pRQhnUjmEepS6UOk\/TKOovu1Hzr3pcW2GIEhi\/8uhGXn1vvsCYYmwabPXMq+FnvxEt1WdVGL1P6Ub9xhtvpJtvvtnqZTZqZM\/GbvGz\/y\/\/8i+eg5L0WY0SSr8XnLAVx0FCaVM+xcR2RbDZt5wIpd8EuOpcYXOUtkKplsQHdUzlcfp1+J4glOZbN4\/x266ONB8qU6TCFmHpYh9XKP2E3Hxw\/MriN7QVJO56f9A9W\/Xwb9u2zfOAx4wZ4w3VBBn2II9SN\/4uhNLPK+U+z2XXvSxlkM2hV7\/Fc37PhJ8oRnnEJoOoOUrmykzV0GuYUIatNOV4XN8LLriAWCz8+nWUMCQVSt2G6SNSprcb5VEmFUq9T82ZM4fmzp0bulXN1g4rzyvq5V0JZdJnlVdT+\/U1v75jvkjqIwbmUGvQ7\/oLlzllZPaRqC1AIkIZtNpKf\/sPMrg2QwBqD6bquEHeWNADoxorzrCm7g2F7aO0Kb8S8rhDr37l9nujDtprpv4e5S34GdMg70VnqJcv7tCrzcOnyqWXRc8nzn4q\/SVLvcXrIxXmm71u2M25Hz\/v0IVQcp7mylO\/+eOwlxddLIPmbIOGn1TepmeuswszUrqhirPq1UYor776avr617\/udQn2dvyGM4M8wzhC6eel6v3v9ttvpxkzZpDfEGmQR2lT1qD+Y74shi3ui2OH49itpM8qv9zxM2buRfcTyrDpjTRCaTNv69fmIkOvuqjoK5b83vTNAkRNKoct7oharGMav6SLefzqxHXjTnveeed1W0Gl6q0aOKqOPIQcNjxqu0hF5avSUvUPWszD4fW3fj8+QUOzulEyBdk0WGET8uZboN\/CBpVe0OpB9X3YEnT9pc1vhWKYYTdF1s97thXKoI33UQvTbJ6loPbXDyngMEHDwGHD2hzPb+45SJzMubigeqs+qB9m4OcZ6C8Rep6utx2EPav6fLZ+qIrZ\/\/h3c84\/aM7Ub35NZ62Xx2aVb1j59fg2NilM0KLmn1VfMfPRF\/PoAhw2p+xCKJkp1\/+b3\/ym96Lj99HLI7KYR2VqM\/4bNtelHiab7SHmIpyopc7m93G2h4QJJU8YR+WtG6egOupceEhQncqSdntI0INoGqiotrMtn\/lgBG0PCXtLjSoLlz1oGN4cfdAfCHNhkZ6PuSiH46k5SrN9v\/\/97xPPPethXAmlKm\/UikP95crsI3q9\/IQkqqx+eZvee1D5op7LKANlir1eN7NfuBZJlbffC4Pf6ly2C2yrFi9eTP\/0T\/9Umj9UL6bm0GZUn9VflvW62Sx+0bnaPD82NilKKNn2ReVlfs\/tqexbWJ9i3mqlsGIfNEcZZkvMtuSwUdtQxLaHBL1Vmn\/3e6v3W0XG8UyDG3efnW2ZKiVcnOGhcpS50stXDiblyFOfswx6dmzKpdozzvC1TboI456AEoiwl2b3ufbcFJVOOT9wwBYpNzh\/eJJX\/zlIUDkcvzmrgjc1NXnnw\/Lv06ZNI57gNoeVbMtSaeEqXYgqvXyV1p6VXh71TOlH2FV6mXti+aL2xPZEJtJ1Vp6\/\/hLpZNWrTcHVg6nEjzsAr+JauHCh74HkZpq6sMaNa1O+coepdCGq9PKVu\/3ymD+3qd+h6HmsSxHLrA9vw\/PPpoXVsCvnJnIoelQ1WNx4jxpvSmcv0Pw9Kr4ulEmHn9rb26OywfcgAAIgAALCBHgeN0+fzDxK0wuMM3yqhh94s7jfClGbK6VYJHm106ZNm\/LUPigrCIAACBSOwIUXXki81SYvglnxQuk3l8LeJa84U66x+btfr1Ljztw4Z555ZuE6XiVUiF9CeOsEGMu0BvjKcNVTBePsGK9fvx5CaeJOMvRqu+DAZhhXCWWeGke+y7rNAYzd8jRTA19Zvpw6GIOxH4HMPEpzqDVqQY650jWs+aLSwgMg3\/nBWJ4xjDgYyxOQzyGP\/TgzoWT8tttD1Moj3syrbl9Xzae+a2xs9LaKhIXVmzyPjSPfZd3mwBuJ77vvPuKrqqqqqtwmjtS8jdrgK9sRwFiWb15fqDMVyrADB\/TD03nhjXmxJgMOOmXG5oQOCKX8A\/DOO+\/Qzp07aciQIRBKAdzgKwDVSBKM5Rnn0RZnKpTyTRCcQx4bp5y8kuQNI5OEmn0c8LVnlTQkGCclZx8vj7YYQmnfvggZQQBGRraLgK8sX04djOUZQyjlGSfOIY+Nk7iyZYoIIyMLHnxl+UIo5flijjIbxolzgVAmRmcdEYbcGlWigOCbCFusSGAcC1eiwHm0xRh6TdTUiORHAEZGtl+AbzhfF0dUHj161FuQdsYZZ2BBWoLubHPSTiGEMury1iB2fIzc6tWrE6DNJkoeGycbMu5ygSF3xxIvIvFY4ojKeLykQtscTZdHW9zNozRv+bABqg4ph1Da0CpuGAilbNuCbzBfHFEp2\/dsUlfH\/0WdfgahhEdp058KGwaGXLZpwTdaKKOMtGwL9ezUbQXQNlwl0QydowzzLpNedVWuyuexccrFKmm+MORJydnFA18IpV1PKU8oWxtrG648tfDPFUJZSa2R87LAkMs2IPhCKGV7WLrUbQXQNly60riN7SuU6jb7qKymTJnS7SzWqDjl+j6PjVMuVknzhSFPSs4uHvhCKO16SnlC2dpY23DlqYVjj7KSKmFTljw2jk29KikMDLlsa4BvsYWSL41oaGjoctnDhAkTvMvq8\/CxtbG24SqpzthHWUmtkfOywJDLNiD4xhPKHW+8I9sgKVMfWntSlxT4JqRbb72VJk2aRK+\/\/jo99thjuRmx44rYCqBtuJR4nUa3Ekq\/odh58+Z5bz55+eSxcfLCVpUThly2xcA3nlDe9tOX6bafbpNtlBSpz7qsgWZd9pEuKfDdurfccgsdf\/zxtHDhQqqtrU2RQ7ZRbW2sbbhsSx+eW6RQskiuWrWKli5dWmq0OJcqV0pl89g4lcLOthww5LakkoUD33hCyR5lJXuV7FGaXiXXUB+CTdZTyhPL1sbahitPLfxzxarXSmqNnJcFhly2AcE3nlDKtoZM6iwijzzyiJc4D8EOHz5cJiOBVG0F0DacQBETJwmhTIwOEU0CMOSyfQJ8iy2UPFJ3880304033ujNUS5btoyuu+46qq6ulu1YjlK3FUDbcI6K5SQZDL06wYhEmAAMuWw\/AN9iC+U3vvEN0le58rTXtm3bcrOgx1YAbcPJPk3xUo8USk4Oi3niQe2poWHIZVsefIstlLK9Rz51WwG0DSdfYvscrITSPrnKDZnHxqlcmv4lgyGXbTHwhVDK9rB0qdvaWNtw6UrjNjaE0i3PHp0aDLls84MvhFK2h6VL3VYAbcOlK43b2IHXbG3ZsiUyJ76DUt82EhmhjAHy2DhlxJUoaxjyRNisI4EvhNK6s5QhoK2NtQ1XhioEZhnpUYbto5w1a1bhjleqpMbJW1lgyGVbDHwhlLI9LF3qtgJoGy5dadzGxvYQtzx7dGow5LLND74QSibAQsOHEkSN5mV9cIGtANqGk32a4qUOoYzHC6FDCMCQy3YP8IVQQihln7Gg1FMNvTY1NeXmvNc8vsWUp0skzxWGPDk7m5jgW2yhDLs9hM+AbW5upo6ODuLrDdmeKY9S2TamM3bsWJo\/fz6tWbOGWlpaPGBZnctta2Ntw9k8E1mFiRRK9QYzceLELmVasWKF8\/lJfb9mWOOqs2bVgiPVOcJOsPjlg\/vou99aRv+w6a+pvr4+K749Kh8YctnmBt9iC2XQ7SHm2dosqEoo+QQfFtAFCxYQL66cPXs21dXVeYcUYOjV3fNoJZTusgtOid+YZs6cSa2trV4g9bN51iF3Ju4MjY2NXe5tU50jKIfHb99Njy\/YTVduroVQCjUoDLkQ2A+SBd94Qnl0b+XeHMI1qRrY0K1CfreH8N\/mzp1buk1E\/\/3555+nlStXel4kOwr6\/OU999xTut9StmceS93WU7QNl0WZbfOoGKFkb3Ljxo2lBo\/zNmTG9as8hNK2SyQPB0OenJ1NTPCNJ5T7fnAT7XvoZhu0ZQnT\/8s3Uv+v3NQtb9P2mYt3dKFct25dF7upfwehdNeskUKpj42b2brcR8mdgz88ZMAf8\/ewKscRyj97+G3Poxw8eLA7ikjJIwBDLtsRwDeY71NPPUU8PbR+\/frSiBF7lEf2VK5X2ev0hm5epd\/tIXnzKH\/+85\/7jtrt2rXLa8D29vZubSX75KRPPVQo9WHOcePGeUOefGjviBEjaPLkyZ6ojR49On0pfO5gsz0Q2PZuTOVR3tvnC155r7zySu8aG3zcETh8+DDt3bvXewmpqqpylzBS8giAb3BH4PUK3\/rWt7oIZd66TdDtIVwPfbqp0ucoH3jgAV9HhG9DWb58ealZ9JeaSm+rWNtD9CEBfvPRx8bTVtQcbrARSiXknLcaow8qhxLK8+\/aT6cO6eU1JLzKtK3WNT63x+7du723SQilW7acGvgGM2V7dNVVV+VaKMNuD9FH9viFgH\/n67hqa2tLc4NMR1\/YqBZHZr3q9Wc\/+1mgR8le5aZNm2jRokW5aqtYQqmLl+2mV1tzEXfoNY5IchmUUF7z1Ll06pDetsVCuBgEMDQYA1aCoOAbLpTm0GsCxIiSgoDtIh3bcCmK4jxq5BylLmC6OJqTyGlLZnqQYYt5lEhGrXTVywShTNtC0fFhyKMZpQkBvhDKNP1HOq6tANqGky5vnPQjhdLcjsECtmTJEm+vTltbG5nbN+Jkroe13R7CcbgMvPE2argVQpm0NZLFgyFPxs02FvhCKG37SjnC2Qqgbbhy1CEoz0ihzLKwQQcOKLHWFxKZt5tErcBVHuVXV59Fwxr7ZlmtHpMXDLlsU4MvhFK2h6VL3VYAbcOlK43b2LHmKPWsXc9Ruq1W99QglNKEsT1EmjCEEkIp3cfSpG8rgLbh0pTFdVwIpWuiPTg9GHLZxgffYgslbw+ZNm0azZkzx1s1qrbjxdmCZ+65lO2RXVO3FUDbcFmWPSovX6HUh0DDEuDDedUBAVEZlft7eJTyLQBDLssYfCGUUT0MQhlFKNn3iT3KZNmVL5YSys\/dMYRGju9fvoIUOGcYctnGBd9iCyXvo1y7dq23UPLuu++me++9lz71qU95t4Tw4kV9P6S+r1Ktz2A6fBAMr99Qf+OzYPULLSQus1CtYusp2oaTfZripV5Ri3niFT1eaAhlPF5JQsOQJ6FmHwd8iy2UfkOvXGNe3c\/iN336dG+nwYABA7qcjKbvAuDj4dQB6nyziH65hM1Rn\/a9sXtIWwG0DZemLK7jQihdE+3B6cGQyzY++MYTyvffbpdtkJSpH39y1+v+\/IRS3ZKkf8cCyOKo7qPUh1v5O\/2mEb2Irk9TM6tvK4C24VLidRrdSij1OUt23bdv397lxHqnJRJKDB6lEFgtWRhyWcbgG08oD239Hh3auki2UVKkXn3ON6j6nL8vpRC2mMcUSvN+YLWvnRPThVLte1eZ2Nzdm7RKtgJoGy5pOSTiRQqlcuvZhZ86daq3eMe8IFSiYK7TVEJ58fRBdPGMQa6TR3q4PUS8D0Ao4wkle5TvVbBXecLJ9aR7lXGEMuicbfOuSt3zhEeZ\/BG1Xsxj3hiS132UWMyTvLNExYQhjyKU7nvwjSeU6WhnH9tWKM05Sh7xW7VqlTcUqw+98kIeJZR8qTNvN+FPnBPN4lCw9RRtw8XJWzpsjxNKeJRyXQqGXI4tpwy+xRZKdQLZ5s2bS6te+TQy3kepiygfG6qvetWPE1XXDjKpu+66i1pbW0sraXl\/5kMPPUQLFy70bh1x\/bEVQNtwrsuXJr3IoVe1UkofelXeZVNTE40fPz5N\/pnFxdCrPGoYclnG4FtsoZTtPfKp2wqgbTj5EtvnECmUnJSqmJ5sVnec2VclPKQSSt5DycOv+LgnAEPunqmeIvhCKGV7WLrUbQXQNly60riNbSWUbrMsT2oQSnnuMOSyjMEXQinbw9KlbiuAtuHSlcZtbAilW549OjUYctnmB18IpWwPS5e6rQDahktXGrexI4VSnzQ2s4662sptUdOlBo8yHT+b2DDkNpSShwFfCGXy3iMf01YAbcPJl9g+h1ChVKuweFVVXg4\/D6q6EsqhjX3oitVn2xNCSGsCMOTWqBIFBN9iC2U5bw+xOd7OXHlrtoatANqGS\/SQCEWy3h4S56oXobKmShZCmQqfVWQYcitMiQOBL4QyqvMkvT0EQhlO1sqjVHt5ohqpkr+HUMq3Dgy5LGPwLbZQZn17iPLseMRwzJgxdODAAe8wAv7w4QR8kwl\/+NhSdRqbut2ED2fn\/Zz6x9ZTtA0n+zTFSz1yjjJvJ\/BEDb2eOqQ3XfPUufEoIbQVARhyK0yJA4FvsYUyy9tD+JaR5uZmWrBgQUkEmS4L5R133OFd6+V3a4m6WNoUSY5rK4C24RI\/KAIRuwmlOtmBr3WJ+uRxMQ+EMqpVk38PQ56cnU1M8I0nlPtfedcGa9nCsC3SP1neHsL2XT8vVp0De8MNN9B3vvMdUqOIap0K32Jy6aWXEoSybN0lm4zV0CuEUo43DLkcW04ZfOMJpXrmZVsleermcZq2Z73yea5pbw9Zs2ZNlxuglFDqJ7CpdSl8XmxDQwOEMqiplXfJK17NxTx5G5LVH5o5u0cm792IGUgAhly2c4BvPKFkj\/LNCvYqa4b0Jt2rjCOUaW8PgUcZ71lNvOoVQhkPdE8IDUMu28rgG08oZVvDfeq2Quni9hAeUp08eXKXaxO5Rpij9G9XX6HUL2oO6w5TpkzJzf5KeJTuH2wzRRhyWcbgW2yhzPr2EP0M7+uvv55efPFFuu666zzI5qpXHlHUy4dVr1pfDBt6lTUJ7lOHULpnCqGUZ6rnAKEstlBm25vc52a7mtU2nPsSJk8xcntI8qQrK6YulLw9xFxxVlmlzWdpYMhl2w18IZSyPSxd6rYCaBsuXWncxs5UKPWtJ7ZbS5S7rx96oP6mNsQykqhhYAil247jlxoMuSxj8IVQyvawdKnbCqBtuHSlcRs7U6HkZcb84VW0+s9BVdIFkU+HUCtvo84c9EsPQum240Ao5XmaOUAoIZTZ9zr7HG0F0Dacfc7yITMTSnO+M+pMQnVryfnnn087duzwxFUJZVRcCKV8x4FQZs8YQhktlPxCXV9fn33jIEfi0354f+f69etD2wBCGdJZWNx4M2tra6t3RqD5uxn11Vdf9f5UXV1dWsashDLJ1hTdoxxz54n0x1\/EDSKun20YctdEu6YHvsF82Ui3tLTQpk2bZBsBqYcSYMeG93j6fXbt2uX92VZQKwl1Yo9SX1rMh+r6LRfWK2p6gbbDp34rb83tKzbznbpQ\/vikb9OlzefRpEmTKqktcl+Ww4cP0969e2nw4MFUVVWV+\/pUWgXAN7xF2BDv3r07VbMx49\/\/\/vfUv39\/9OEEJAcNGuQ9\/36fZcuW0fLly0tfRXmeCbIXi5JYKOOWyKVQ8vymOrSXPU7zd7+ymR5l3R+fGNigceuG8McI8JwyGyoe+oJQuu8V4OueqZkiGMsx5hcZ\/sde\/6JFiyKHaOVKEj\/lTIUyztCrqorNXs6oYVxOi4XyJwufpX6dp9Pn7hhCI8f3j08LMUIJYGhQtoOAryxfTh2M5RljjjKEsTnUarsgx1Yo586dSwsXLqTa2lrfUkAo5R8AGBlZxuAryxdCKc+XcyikUKrVpzzUaX5s5gb1OHG3h3BcUyj1a1\/Gjx9fOlaJ50l5ZWzQB0Ip\/xDAkMsyBl9ZvhBKeb6FFEpdlMaNG+ed\/8cb\/0eMGNFtJaoN4rADB9Q1L3woL887hg29mgcOjB071jvMV49nlkcJJf99wqLzMfRq02Axw8CQxwQWMzj4xgSWIDgYJ4AWM0rhPErTm1P3krEnFyRsMZllFlwXyr+a9km6eMagzPLuKRnByMi2NPjK8oVHKc+3kB6lKZS8LWPbtm3eEGeSvYzZNIN\/LhBKefow5LKMwVeWL4RSnm8hhZIrpc8r6uK4bt26LjdkZ4M4eS66UF70lXO8la\/4uCUAQ+6Wp5ka+MryhVDK8y2sUJqLZ1g4lyxZQjaHDGSD3S4XJZR8q\/jHP\/0RCKUdtlihYMhj4YodGHxjI4sdAYxjI4sdoXBzlLEJVHAECKV848DIyDIGX1m+8Cjl+RbSowzbw5jXOUp4lHIPAwy5HFsYcVm2KnX0YXnOhfMoiyiUwy7qQ33fP52uWI1D0V0\/EjAyrol2TQ98ZfniZUSeb6E8SvPQ8SB8UZclZ4PdLhc19AqhtOOVJBQMeRJq9nHA155V0pBgnJScfbwe5VHaY6mMkBBK+XaAkZFlDL6yfOFRyvMtlEeZDa5sc1FCObKpP72xoQ9d89S52RagB+QGQy7byOAryxdCKc+30EKpD8XyDeLbt2\/P1R5KbhwIpfxDAEMuyxh8ZflCKOX5FlYo1V2PfEXW1KlTvVN5+DB0Pvc16iDybLDb5aIL5csre9Gc3SPtIiKUNQEYcmtUiQKCbyJssSKBcSxciQIXeo7SPAg9r9tD+IzXZ1o7IZSJunh4JBgZAahakuAryxcepTzfQnqU+vYQCGU2nSjPucCQy7Ye+MryhVDK8y2kUHKleH5y48aNpA+9KtFsamoivkkkDx819MoeJQ+9YjGP+1aXnhSdAAAgAElEQVSDIXfPVE8RfGX5Qijl+RZWKPWK6RjnzZuXG5HkciuhHHvHEHr0fx32hPLUIb2z6Rk9JBcYctmGBl9ZvhBKeb6FFsps8MnmAqGU5QsjA77yBORzwMuIPOPCLeaRR5ZdDhBKedYwMrKMwVeWL1725PkW1qN84YUXqLm5mTo6OrpR5G0iS5cupdra2mwIp8hFCeUVq8+ih79wkL66+iwa1tg3RYqIahKAIZftE+AryxdCKc+3kEKp7qLM037JoKaGUMo\/BDDksozBV5YvhFKebyGFMuz2kGyQustFCSUv4ll+\/hvexc0jx\/d3lwFSIhhy2U4AvrJ8IZTyfAsplMqjnDBhAo0ePTobikK5QCiFwGrJwpDLMgZfWb4QSnm+hRRKVSk+xi4vc5FRQ6\/wKOUeBhhyObYw4rJsVerow\/KcC7HqVQ23btmyJZJYHhfzsFDyYh4eduXDB\/BxRwBGxh1Lv5TAV5YvXkbk+RbWo8wGnXwu+tArCyWLJOYo3XKHIXfL00wNfGX5Qijl+UIos2GcOBdTKIc19vEW9ODjjgAMuTuW8ChlWQaljj4sz70QQ686prBVrxK3h+j3Xtoekcfzpw0NDZHH6elCyUfY1QzpDaF0\/EzAyDgGaiQHvrJ84VHK8y2kR5mlUPLBBnzwemtrq9da6ufhw4cHth6L5JIlS8hGVCGU8g8BDLksY\/CV5QuhlOdbKKHUPbswdFOmTPEucnbxUbeUzJ8\/n6qrqynMU1QCXlNT42X92c9+NpZHyfdR8gdDry5a7sM0YMjd8jRTA19ZvhBKeb6FEkqFK8sDB1gY+aOE1\/xdb0Iu1759+4hPDJo9ezY1NjZaC+Xnf9SXXltxJr328tt0xeqzs+kZPSQXGHLZhgZfWb4QSlm+u3bt8jJob2+niRMn0vr166m+vl42U0epH9fZ2XnMvSrzx\/Qg2cPctm1bqMeqDkSII5QPnvx3dMnha+mikZfSn\/\/jyWWudbGyP3z4MO3du5cGDx5MVVVVxapcBdQGfOUbAYzlGC9btoyWL19eygBCmYB1VkJ5weL99O7jDfTGhj40eUNDgpIiShABfnHZvXu395YIoXTfT8DXPVMzRTCWY8weJf\/btGkTLVq0CB5lEtRxhl5V+kk8Sj5w4OWVveiXD+7zLm\/Gxx0BDA26Y+mXEvjK8uXUwVieceG2h8gj+zAHc6jVZttHUqF844k+9Mi1r9Cc3SOzrGLh84KRkW1i8JXlC6GU58s5QChTcE6yPQRCmQK4QFQYcgGoWpLgK8sXQinPt8cJpXor4Irz6tO2tjYK2\/No0wRBBw4E3WICobShml0YGHJZ1uAryxdCKc+3xwllNkjd5aIfOPDejhp64PKXvDnKU4f0dpdJD08Jhly2A4CvLF8IpTxfCGU2jBPnAqFMjM46Igy5NapEAcE3EbZYkcA4Fq5EgQs3R5nlEXaJiMeIpAtlv87T6f7LX6LP3VFPwxr7xkgFQcMIwMjI9g\/wleULj1KebyE9yiIL5dprX\/Gu2oJQuns4YMjdsfRLCXxl+UIo5fkWSijLcdardBPpHiVviJ876JfeWa+4k9IdeRhydywhlLIsg1JHH5bn3qOGXuVxus0BQumWJwy5PE8zBxhxeeZgLM+4cEIpjyy7HEyhXHzBb72hV3iU7toARsYdS7yIyLKER1kevoUaejUR+g3F2twBWb6m6J6zn1CySLJY4uOGAITSDUcYcVmOYamjD8uzL6RHySK5atUqWrp0KdXW1noU1SKfpqamyOut5LHb5eAnlMMa++BOSjt8VqFgZKwwJQ4EvonRWUcEY2tUiQMWTiiLuuqVF\/Pcf\/mLVDOkN4QycXfvHhFGxiFMn6TAV5Yvpw7G8owhlPKME+dgepQQysQoAyPCyLhnqqcIvrJ8IZTyfDmHwgklV6qoQ698e8ibr7xLV6w+O5ve0QNygSGXbWTwleULoZTnW1ihVGLZ0tLShWLeF\/NAKN0\/FDDk7pnCo5RlaqaOPizPu5AepTy2bHIwh175d1ze7JY9jIxbnjDisjz9UkcflmcOoZRnnDgHP6F8fMFuXN6cmGj3iDAyDmH6JAW+snwx9CrPtzBDr2ql65YtWyKpjRo1qsu2kcgIZQxgCiV7kzz8Omf3yDKWqlhZw5DLtif4yvKFUMrzLYxQmqjCFvPMmjWLRo8enQ3dlLlAKFMCtIgOQ24BKUUQ8E0BzzIqGFuCShGscEOvRd5HuX3jQVzenKKz+0WFkXEM1EgOfGX5wqOU51tIjxJCmU3HKUouMOSyLQm+snwhlPJ8CymUXKmi7qPc\/8q7xAejf3X1WbiT0tHzAUPuCGRAMuAryxdCKc+3sEKpV0zHuGLFitzMT3K5zTlK\/hvupHT7YMCQu+Vppga+snwhlPJ8CymUhw4dIv6nDkPPBqNMLhBKGa56qjDksozBV5YvhFKebyGFMuqWkGXLltHYsWNzIaR+Qok7Kd0+GDDkbnnCo5Tl6Zc6+rA888KtetXVnwVx\/vz5VF1dTS+88AI1NzfTwIEDc7uPkuvGQok7Kd09GDAy7ljCiMuyDEodfVieeyGFkrEpz3Lv3r2eOPJhBFOmTCHeR5mXT5BHiTsp3bUgjIw7lhBKWZYQyvLwLeTQq45SeZEdHR2UtwPRuR5+Qomrttw+LBBKtzwx9CrLEy8j2fMttFDedttttGTJEs+LvOSSS2jixIne3KQairXFrR+PZ3P8ncqX09dX2fICo9mzZ9PatWtLWUd5uH5CyUfY8edzdwyxrQLChRCAUMp2D\/CV5cupg7E848INvepDrm1tbTR8+HCPovo7\/7x06VLrxTwsfPzhIVv9Z7+mYZgchtN\/\/vnnSz\/zClzOf9q0aTRnzpxSmaKaN0gocSdlFDn772Fk7FklCQm+SajFiwPG8XglCV1IoWSvbdKkSb484qx6NU\/54aHcuXPn0sKFC32FVhdS5UFOmDDB27sZFdevsH5CyQej89+veercJO2NOAYBGBnZLgG+snzhUcrz5RwKJ5QusbG4zZw5k1pbWz0v0Pxdz0sJY2NjI40fP97by8lDrep33du03eOphPLzP+pL9fX1NHjwYHruhwe8G0Rmvvoxl1XtsWnBkMs2PfjK8oVQyvLdtWuXl0F7e7s3fbd+\/XrPFufhc1xnZ2dnFgU1vcCw4VPTg+TysYfZ0NDgCScfq9fS0lIqts18pxLKB0\/+Oy\/elVdeSX962kR68h\/eoQlP9ssCQeHzOHz4MPHKaH4JqaqqKnx9s64g+MoTB2M5xjwCuXz58lIGEEof1i6FkkWTV9+qxUTm735NrYTygsX76cwzz\/SM+eEX+9IPJ3bQVRsaqGZIb7ke0kNS5hec3bt3e2+JEEr3jQ6+7pmaKYKxHGP2KPnfpk2baNGiRfAo\/VC7HHo10w8bxlVh\/eYo+aqtR65tpytWn0WnQihTPyEYGkyNMDQB8JXly6mDsTxjzFGGMDaHWqMW5OhDrX5DsXpWUWlxWD+hxA0ibh8KGBm3PM3UwFeWL4RSni\/nAKGM4Oxiewgfoacv7FEiWldXF3pSkJ9QcnFxg4i7hwOG3B1Lv5TAV5YvhFKeb48TSvVWwBVnkdL3WQbhDjtwgNNbuXJll0MMbA8csDn8AEIp\/xDAkMsyBl9ZvhBKeb49TiizQeoulyChxA0i7hjDkLtjCY9SlmVQ6ujD8twLN\/Ra9PsouUvgvFd3DwaMjDuWEEpZlhDK8vAtpEdZ9PsoudHY0+RFPTjvNf2DA6FMzzAsBfCV5YuhV3m+hRRKvVJFvI9SCSVvE7li9dnZ9JIC5wJDLtu44CvLF0Ipz7ewQskVK+p9lEoo+cxXnPea\/iGBIU\/PEB6lLMOo1NGHowil\/75wc5Q6kiLeR8n1Y5Hk817n7B6Zvgf08BRgZGQ7APjK8oVHKc+30B6lq\/sos2kG\/1yCVr0qoWSPEqfzpGshGPJ0\/KJig28UofTfg3F6hlEpFM6jdH0fZRRAye+DhJLnJx+4\/CVv6BVCma4FYGTS8YuKDb5RhNJ\/D8bpGUalUEihdHUfZRQ86e+DhBLH2LkjDyPjjqVfSuAryxdDr\/J8CzP0al6wbINO3Q+5evVqm+BlCRMklFwYHGPnpklgyN1wDEoFfGX5Qijl+UIob7uN8iqUfDrPyPH96eIZg7LpKQXNBYZctmHBV5YvhFKeb+GEcsuWLbGo8eXJeRVKPp1nWGNfCGWsFu8eGIY8JcCI6OAryxdCKc+3MEKZDarscwkbeuXtIfzB6Tzp2gWGPB2\/qNjgG0Uo\/fdgnJ5hVAqFW8wTVeE8fR8llG++8i5O50nZoDAyKQHCo5QFaJE6+rAFpJRBIJQpAUpGDxNKHDrghjyMjBuOQamAryxfDL3K88XQazaME+cSJpRqLyVO50mM14sIQ56OX1Rs8I0ilP57ME7PMCoFeJRRhMr4vY1Q4tCBdA0EI5OOX1Rs8I0ilP57ME7PMCqFQgqlvq9yxIgRNHnyZOIVsbzKdenSpVRbWxvFpSK+DxNKLiD2UqZvJhiZ9AzDUgBfWb4YFZHnW9ihVz7nlT+zZs2iBx98kFatWuUJ5Lp162jbtm3e3\/PwiRJK3kv5uTvqvW0i+CQjAEOejJttLPC1JZU8HBgnZ2cbs3Aepe5Nsgc5e\/Zsqqur88RRncaTF6\/SRij5wAE+eACfZARgZJJxs40FvrakkocD4+TsbGMWWijVsGtTUxONHz++cEKZZC\/l0b3b6MiebaX+0ev0Bqoa2GDbXwoXDkZGtknBV5Yvhl7l+RZy6PXQoUOeF9nY2EjDhg2j6dOnU1tbGw0fPpz0Idls8KbLJcqj5O9tL3A+9OtHiUVy3w9upqrTG+joB2JZ\/fExdGTvNqr+2CXU\/ys3pStwDmPDkMs2GvjK8oVQyvMtpFBypfQLm6dMmeINu7JIdnR00Pz586m6ujobuilziRJK272U+35wEx3du90rTf+v3FjyIJV3qQS035hJH4TpOYIJQ56yk0ZEB19ZvhBKeb6FFcps0MnnEiWU6rqtsC0iBx69jw48uizSY2SxZK+TPzbh5WufTQ4w5LKcwVeWL4RSni+EMhvGiXOJEkpOmLeIfHX1Wb4rX1n4dt70GRr6jy9bz0MqL5PjDbymjXhotshzmDDkibunVUTwtcKUKhAYp8JnFblwi3lUrXlbSEtLi\/frihUraPv27bRx48ZCDb1y3YKu22LB27O4mXg4td+Yv7XqDGYgFtp9D92cKo1EGWcYCUZGFjb4yvKFRynPt7AepZqPnDlzJk2dOtWbozS3imSDt2suunjPmzfPW4kb9rHxKHnl63uH2umv5h1HvU4bXUqO5yUPPfcY1d30i1RV5aHbd379mJdG1cBhhVvwA0OeqntERgbfSESpA4BxaoSRCRTOo\/Q7lYeFcvTo0WXdHsILjFi4W1tbvUZRP\/Nq3KCPEsr\/vuEC+pOP1nQJ9v7b7XTk9SfpmbZn6T8Wf4m+uf7n9P7br1LvIV+id1\/+He3\/8VI6\/Zq2VMOmPAfKN5TseOItb8Use5h79v+RJ5j8UQcdDL2oT5ffI3tdBQWAkZFtDPCV5QuPUp5vIT3KShVK9ib1oV\/2ehsaGkK9ShbKR+7eSc9N66S11\/xhqUewSB5+5f94wvjcv4+gdQsuJF7QcxKtoUNbFxF\/T+8Pptov\/GeiXsQCuX3jW7R\/B\/9\/kGqG9KahF\/X1fu9zUocnmByGxXPU926nZ1ft65LPqUN6e\/FGNh07KvDUIb0q9vQgGPJEXcQ6Evhao0ocEIwTo7OOWDiPkmuuREkfejUPH7Am5CiguYfTZk+n8ih\/df0A+tO69+ibYz\/llea9nT+it156kKrP+YY33Np65nM0\/qGhnhgd2bqBXl9zMw344ufpyGub6KRPzvPiHH9yvVVNnvvhAdrxxEE6evQofbLp2Ik\/+hF5Hb8\/Si9t\/S09v+Hf6C\/+bwttGfJ5+pOz+9Piw5d4YVdv\/igNrT2JPtH7ZOr87XvU+Zv36JVBR+iiI\/3ouD84gU45eDzV\/OXJXpi+522hE\/ZeQP0G7PH+hX1O8Cn\/+70H0\/Hv7vrwBeKto8FJ9K8n2tfe7fujR47Szl076YzBZ1BVr6oPv+fw+KQmACOeGmFkAmAciShxgF27jtmX9vZ2mjhxIq1fv57q6\/NhG47r7OzsjKq5egPQw9nMC0alm\/R704NkMY86d3brPT+mnyx8lm4\/aze9c+44uvaUx6h5wmVU3d5K7\/ceRIfqj51Zu+byt+j0806g0defRJ3fn0BUW0\/HfeV2qnrrGTpx93IvzJH+l9HRPqOIxcXv89bO9+nZe9+lPU+\/56XD6fGHhXHngaO0ZNOb3v9n9Kvy\/ucP\/zx28JueAH1h97Ey7XlnNx3ft4p2Vg2m4\/ucQB1Vg+nc\/XupblgDPb2xhk49ePyxdJ\/\/JB147XTv57pzfkVbN\/6Z9z9\/WDQPvH566fv9fd+jc059lt46XOdb9mNCu9tLkz\/7t9fSqcPe8OJzOio9Jcb8N5UP\/6\/C6D\/rf+t32h4vLVVOVYhzGn\/ulVulr5dd\/5uqY8fWT5TKz2mp8iXtU2Y8VT+9vHr9dd56WVR4lZ7+O\/+swnqcP2Chs1LtVgp32rEXH9WGXl015hyO01Ef\/l3xUGFVXfSw3s8D9pTaQf3MbX\/g9UFd8tPLoufDbdYt3w\/6mqobl5XLYfYXr86vne6VvW7Es6U+rBiptP046HFVOygGeruofqzy0cPq4fR2VfVRcXSeqp0uGLvSK69iqPJWv+t9VNVHtYnXr402NfuI90x\/0I5mXYP6l\/r7OY3\/0aVsKj8zPfV3PzYmA5Mb89LbXvVNVVfz2VF8uJ25b3nP+bk\/p8sf+OtiCaUr4+MqnSRCyQtyfvPArbTzylX06JuDqO6UKpr6sQ2ep3jyhf9cKtp\/fvd1byj0kslb6dDa26h29r92KfaR156k9357vSeUJ57SQO\/3GkwnDv3rUphfrnqD2JM87SMn04jZ\/WjHG4foiRffpCe3HaSqqmNe1pfOPex5gef3+ZU3N8oeXtWA0XT09Sfp+Hd3U+f+p6jXgNH03tvt3ncsmu+9uonefvJX3lBwnwuuoP0\/vtfL+6T6vvT+waPH5j2P\/2M63K\/R+7nfR3vTi0+N9E4O4k+vD47We+udOm\/Il\/\/3M+bnfPo578\/8PXvOPESsPqcO7f1hnA+8Vt1o723vRwcPHqQhH32XDu0\/0wvrJ3yeIfnAOJbK8IEY2PQRUzjNfGzSiAqjvxDoddR\/1l8ozPB+LxS6WJgvJKoOx\/7\/UKz0OL36bKeDBw5STf8aqqrqFRhOH1FQ5TXTZ0PvVx42ZJw\/f9TLkhmOv\/MEcMSzpXAl4\/jBS5YyiuaLkx5H9QPO88N+cOwlLYiV4qbE3G\/0ROVZd86ztHXjn38wyvIhU71sZjmZ8fObG6impj8N\/fhvupRFr7eZv15PFU6vg8pH73cqjtmP9PjqOz2+X5v65RnUvn4vs355mn00rE1VG\/q9TKs6qPL8+49fop8+ejd9b+0\/QiijDFGa75MMvSqhHLr4Zdr2Tl+a96PH6X9f9GtPsKrP+ftScdQJPV+74bvecXRBK10Pbf2eN6\/JQseixgt1Xto6x5tPvHD8BqI\/OIGe+N2bdOZJr1FD\/UdpW\/vv6MKaraV8ep12oSfSahi06rQLvXTUJ2h4l0Uw7FPOfZoYtkrTq6Pjgm80o7QhwDgtwej4hZuj1Bfz8EpX\/VPO20PMoVabxTy6UPK4+LzVG+jr\/Vrol4Pupj\/95IfDeFxH3k95fv+\/oQtWPBS50pU9PBbJ1\/\/rcbp3zij6TPMib6ip\/qTXSnOZx1fXEwsjCysLIv\/OAmk71xnd9SojBIyMbDuAryxfTh2M5RlDKOUZezkk2R5iCiV7hI9u+jdq2XEjPfPtT3cp+f2Xv0gfGbSGLvr+NyNrxMO0z275PT0wZzv96qzD3jVd7KXOuuwjkXGLFgBGRrZFwVeWL4RSni\/nUBih1Dfzh6FTh6Rng7drLnEPHNCFsq6W6PcbJ1CfT91Of\/T9Trro7BpaPOEPShk88T++Sy\/8ahT97YY\/C60ai+SK5e3erSMHruhF\/+OLQ7vt0SwHm3LlCUMuSx58ZflCKOX5FkooFa6woddskLrLxRTKg\/81w9sSsunNc+mJF4\/tXWQvkPc17l3cTPf\/8GH63B1DfC9y3vHGO\/TjR\/d4ex4\/ceLJdPH0QfTJUae4K2xOU4Ihl2048JXlC6GU51tIocwGWza56EI56KR2euuZGVTzF497mW\/43Zt0209f9n5ePvBh78aPX1at935nsVQfFsiV\/2+nF77xB709ceT9kfreyGxqU5m5wJDLtgv4yvKFUMrzLaxQ6vdRmhj5zNelS5dSbe2xU2Mq+aML5SnbZ3iLafr84e1dRPCff\/Ik\/c3GZnqkz2V03Kfn0Lbpb9DFMwZRzaXVNHXVb7wh2v2vHKGPv3Qi\/cnZNTSyqT\/xyTn4HCMAQy7bE8BXli\/6sDzfQgrloUOHaPbs2dTY2Ejjxo3zfp4wYQKpk3nUua\/Z4E2Xiy6UfX87wRt2PXHIh\/sfOXU17PrS36ymhU93evsp\/9t\/9iXeP\/hErwM04Y\/OoDd\/9rbnQUIku7cHDHm6PhoVG3yjCKX\/HozTM4xKoTCLeVRFzTlKfRsGV3blypW5uWpLCWX97Y8TC2W\/xpVdbgjhOvPtHvt+cLN356T+4VWw7EnyOat85iqvbMUHQpl1H4ARlycOxvKMCy+U+v7Fcu6jTNKUulCesuPD+UmVFm\/k77jxwwuWk+TR0+PAyMj2APCV5YuhV3m+hRx65Urpp+Do4rhu3bpcXd6shHLgtX9PAzufLC3kUV2Dh10PPrqM+o6ZRNUfH5NNjylYLjDksg0KvrJ8IZTyfAsrlPo8JV+OzMK5ZMkSqquro7a2Ngq7AzIb7Ha5KKHs\/9XR1DDqS12OreMUOm76jHce6sBr2uwSRKhuBGDIZTsF+MryhVDK8y2sUGaDTj4XJZTnTj622lVfyMPDrnsWN1P1xy6h\/l+5Sb4wBc0Bhly2YcFXli+EUp5vIYWyUs96TdKcLJT7\/\/VWOuXyem\/YVT9nlYddd970GTrjpl9g2DUJ3A\/iwJCngGcRFXwtIKUMAsYpAVpEL\/xiHp1BHhfzHFh\/G\/X9y8FUO67rqlYW0aN7t2PY1aKThwWBkUkJMCI6+MryhUcpz7dQHmUeznqN26Qshoe2LqKTP\/2Jbgt5eH6y35hJ1G\/M38ZNFuE1AjDkst0BfGX5Qijl+RZKKBWuop316ieU6pABXsSD1a7pHhQY8nT8omKDbxSh9N+DcXqGUSkUbug1qsJ5+p49yiP7l9BJHxtDp1y0slR0z9N87rHAC5rzVMdylxVGRrYFwFeWLzxKeb6F8yh5+PXOO+8sbQFR3uWWLVs8mvPmzSPeLpKXDwvie0eXUp\/zp3XZGsI3hWDvpJtWhCF3wzEoFfCV5QuhlOdbKKFk13j69OklkVR7KXnvJJ\/vqkSzqakpN2LJQtl50rIuW0N4W8iO\/\/kR78i6qoEN2fSSAucCQy7buOAryxdCKc+3MEJpHjDAFVM3iCxYsIBGjx7t0WSPc+PGjbk563X\/T27xPEr9jFc+25Wv1Kq76RfZ9JCC5wJDLtvA4CvLF0Ipz7cwQum3gMf0MFVl+ZSevFyzdWD9fDry5pIueyh52LVq4DAcMuDo+YAhdwQyIBnwleULoZTnW2ih9PMe87aP8q3NC+jwq4vp7fMep\/r6+tKVWnU3\/wLDro6eDxhyRyAhlLIgQ1JHH5ZHX4hVr2role+d5GFWc35SYWRvsqOjIzdDr+xRbn\/pmEfJQslzljzsal6pJd9NipsDjIxs24KvLF94lPJ8C+NRmvOP7e3t1NzcTPr8pHojWLFiRWnOMhvEyXPRL24efOJRbxEPjqxLztMvJgy5W55mauAryxdCKc+3UELJlVG3hPDPaiuI8i43b96cq5tDuA66UA7Y\/zt659ePYm7S8XMBQ+4YqJEc+MryhVDK8y2cUGaDLLtcdKF8f+7F3pF1uCnELX8Ycrc84VHK8sSoSPZ8IZTlYW6dKwvlvodupv0f\/QsadOJRbAmxJmcfEEJpzypJSPBNQi1eHDCOxytJ6EIs5klS8TzE2bD2QXr2u1fTmKtm0\/DLvoqVrgKNBiMjAFVLEnxl+WLoVZ4vPMpsGCfOJY9vMYkrW6aIMOSy4MFXli+EUp4vhNKCsX5e7KhRoyIPK9AXFOkrbNWiorVr15ZynTJline8XtAHQmnRQCmDwJCnBBgRHXxl+UIo5flCKC0Ys\/DxhwVN\/9kvqn6gwfPPP++FV6cAseBOmzaN5syZQ8OHD7fImQhCaYUpVaBt27bRfffdR1dffbW3VxUftwTA1y1Pv9TAWJ5xHm3xcZ2dnZ3yaKh0kDqLJB9kwOfHzp07lxYuXEi1tbXdiqALqXkIQlTcIOGdOHEirV+\/HkZcqMHz+AAIoRBJFnxFsHZJFIzB2I9AZkLJ4jZz5kxqbW31vEDzd71w5sHs5u9Jjs\/TD0mAtyPzMPDhFPwywsPkYOyeMfi6Z2qmCMbZMc6T05KpUOoeZNjwqelBctOxh9nQ0OBd68Vnz7a0tJRa1Ga+kx+AGTNm0KZNm+R7AnIAARAAARAIJHDhhRfSypUrc0Mol0JpnjNre+4siyX\/wwcEQAAEQKB8BHjEKU+jTiJCqa9u5abgI\/DOO+88Z0OvZvOGDeOWrysgZxAAARAAgSIQEBFKPzDmUGvUghx9qNVvKFbPIyqtIjQU6gACIAACIFAeApkJJVfPxfaQ6upqmj17NjU2NnrzlUHXgJUHJ3IFARAAARAoGoFMhTLswAFelcqTu\/PnzycWQyWsS5Ys8X4OO3Bg7NixubkXs2gdCPUBARAAgaITyFQoiw4T9QMBEAABECgeAQhl8doUNQIBEAABEHBIAH5PsDYAAAUNSURBVELpECaSAgEQAAEQKB6BHiGUQYerF685y1MjczsQl4K3BPFiK3zSE9BXgKvU4l4wkL4UxU7Bj7E6zUvVvK6ujtra2qzPly42MbvambbBXE+Sl35ceKEMO1zdrqkRKooAtudEEUr+vXrJM1884qwgT557z4gZxJhPAOND0sNuJeoZhJLVMugoUn7hUEzz0o8LL5Rhh6sna37EMgn4rVgGpXQE1Jt2TU2Nl9BnP\/vZkoeuvrO9YCBdSYobO4wx19rPyywujWxqxi8fGzdu9HYpsJBOnjzZE02bizKyKaF\/LoUWyqjD1csJvkh5651fbe0pUv3KURc24vv27SN++9b3DXNZ4lwwUI6y5yXPMMam7chLnSq9nLqt4ONEbS\/KKHe9eoRQTpgwwXtjwVuiTHfT54A5B+xrdcfZz2CbQ91J7md1V8L8p+THGPPu7ttVMW1qavJGR\/LUjyGU7vtDj0rRPBkJJyW5bX4IpVuefqkFMW5ubqYFCxaUhgX13+VLVawcFGOulTpUBkJZIW2ModfyNESS+0LLU9LKzzXIiOdlyKryCVPpGEx1LGZQmfX1DnmoV6WU0U8k8zaFUGiP0hxqjTpcvVI6Vt7LgcU97lowaFhw2rRpNGfOnNIl6Ppdr+5y7xkp2c5HYnFP\/P4QNsIU96KM+Lm7i1F4ocT2EHedxS8lcwWmOQ8hm3vxUw8y4nlZVp+HFvJjbI6K8O\/Tp0\/HPsqYDRp1V3Be+nHhhVJ5lX6Hq8dscwQPIGAufJgyZQr2njnqLUFCmZeN2o4wiCYTxNg8cEC\/mEG0QAVJ3G9BFFdt1KhRtHTpUqqtraW89OMeIZQF6XeoBgiAAAiAQBkIQCjLAB1ZggAIgAAI5IcAhDI\/bYWSggAIgAAIlIEAhLIM0JElCIAACIBAfghAKPPTVigpCIAACIBAGQhAKMsAHVkWjwCfMsInt3R0dARW7vbbb6f777+fWltbM72qyXZfK26BKV6\/RI3cEIBQuuGIVECgC4FKOcUl7jmwuFoKHRkEuhOAUKJXgIAAgUoRyrjCF1dYBdAhSRCoOAIQyoprEhSoCAT8hFK\/HmvAgAHeXXzjxo3zNl\/zkC1fqdXW1kZPP\/00tbS0eBjMm1jMTfDmhc46O\/PUJP5Oba5fu3atF1TlOXz48FLUuOJahPZCHUAgjACEEv0DBAQI2AolZ61OKVHXlSnxM48DZAFbtWpVt1NN1LVFZjX8Dqc3y+V3lygOtRfoEEgy1wQglLluPhS+UgnYCqUucqZA6UersefJFzjrd6ty3cNEzRRB2yvQzIuhK5UxygUCWRGAUGZFGvn0KAK2Qjlr1qzSpeJhQnnppZd6Q7VbtmzpxlE\/O1P\/MshbnDhxYuCwK3\/hN2TboxoPlQUBgwCEEl0CBAQIuBbK8847z9t+oi4Stimyn1CqePpcpTlPCaG0oYswPYkAhLIntTbqmhkB10KpPMqg+Ui\/ioUJpSmY+qXFGHrNrJsgo5wQgFDmpKFQzHwRcC2U48ePJxa+O++8s8udiGH3\/ZkHCPjNUfqJIhbz5KuvobTyBCCU8oyRQw8kICGUjJHFUm0d4d\/N7SM6aiWM+gIgvzsCzXsWK2UPaA\/sNqhyhRKAUFZow6BYIOCCQNw9kThwwAV1pFE0AhDKorUo6gMCGoG4whdXWAEbBHoCAQhlT2hl1LFHE8Ch6D26+VF5BwQglA4gIgkQAAEQAIHiEvj\/7gZkHpqpZogAAAAASUVORK5CYII=","height":275,"width":458}}
%---
%[output:97e24fca]
%   data: {"dataType":"textualVariable","outputData":{"name":"xb_max","value":"0.0050"}}
%---
%[output:0eb22992]
%   data: {"dataType":"textualVariable","outputData":{"name":"theta_max","value":"0.0044"}}
%---
%[output:35a48901]
%   data: {"dataType":"textualVariable","outputData":{"name":"disp_percent","value":"5.3008"}}
%---
%[output:534cccc8]
%   data: {"dataType":"textualVariable","outputData":{"name":"rot_percent","value":"0.2820"}}
%---
%[output:4aab1181]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"Ke","rows":4,"type":"double","value":[["-0.6285"],["-4.0121"],["12.0597"],["72.7171"]]}}
%---
%[output:4ac57f12]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_o","rows":4,"type":"double","value":[["-0.0020"],["0"],["0.0873"],["0"]]}}
%---
%[output:724f0e11]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"3"}}
%---
%[output:1fafbf59]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":3,"type":"complex","value":[["-0.6994 + 3.5271i"],["-0.6994 - 3.5271i"],["-0.3800 + 0.0000i"]]}}
%---
%[output:6c24be73]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"A_hat","rows":3,"type":"double","value":[["0","1.0000","-0.0331"],["-11.8839","0","0.3403"],["47.6759","0","-1.7788"]]}}
%---
%[output:720a56da]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B_hat","rows":3,"type":"double","value":[["-0.3992"],["-2.3768"],["37.1102"]]}}
%---
%[output:407f6c67]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"F_hat","rows":3,"type":"double","value":[["0"],["-0.8402"],["4.2977"]]}}
%---
%[output:26003d81]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0_omin","rows":3,"type":"double","value":[["-0.0020"],["0"],["0"]]}}
%---
%[output:59de8216]
%   data: {"dataType":"textualVariable","outputData":{"name":"Ts","value":"1.0000e-04"}}
%---
%[output:2a3ae915]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Integer operands are required for colon operator when used as index."}}
%---
%[output:8f59b930]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Integer operands are required for colon operator when used as index."}}
%---
%[output:8e97ee7e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAcoAAAETCAYAAAC2rUY7AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnV9oHdedx38Ogo2WbUmVBoQsJ1628rJPojFLVbUlpSWkL2KhD5WlB6lCCYYSlmXd6I+bUvxQ60\/sB9VZFmGEI7E40UspcQm4faoahPXgELVPjcqu41VVCSXK0pYogmW1nFHP7bmjmXtnfvfoSPfOR2AS2TO\/mfuZ7\/197jkzd+bU\/v7+vvADAQhAAAIQgEAigVOIkmRAAAIQgAAE0gkgStIBAQhAAAIQqEAAURIPCEAAAhCAAKIMn4F79+5Jf3+\/tLW1ya1bt6SjoyP8TrBFCEAAAhComQAjypoRHi6wu7srP\/rRj+SFF16QDz\/8UO7evSsvvvjiEWyJkqEILC4uyvj4uHR2dsrc3Jy0tLREmzbHemxsTO7cuSO3b9+Wrq6uqru0s7Mjw8PD0XJuraor5ljA7m\/SKvY1\/PznP49eU6Vl3nvvvegDX9qP+SD49ttv59gzFoVA\/RFAlBWO2dramly9elWuX79eaoy2ya2urh5qmkmlzMjy\/fffl97e3vpLB3tcIuCKZ2JionQ8ESWi5G3S+AQQZcoxNpIcGhqSJ554ouxT\/9TUVLTG6OiouP+fVMY01xs3bjD12gDvI1eU7nS6RpQhcNj9daUe326WZdx1QoyEQ7BhGxDISwBRJhCzDeTixYtiRoR2esw2CiNJM8XmjjjtFFX8nKQ7DWun6\/IeJJY\/fgI2E1\/5ylfkl7\/8pZhsmBykidJ+0NrY2Ih23p2yTRKOW8e+WrsN+3t8OrXSVG8WCWZZBlEef\/bYg+MnUEhRxqdDze+vv\/66TE5OSnNzs5hp1XPnzkX\/NaNGK0rT\/EZGRmR6ejq6OCf+uz2cphFeuXJFfvCDH0R\/Zf8fUR5\/4LV7YKXy\/e9\/X959993SOUkjwPg5SnshV3xb9kPU448\/XnaO0mTO1oivY0eEaecc02RZ6RylXSfLMohSmxjWayQChRSl\/fT+8OFD+frXvx4Jz0rSPbim4cVF6Z6zNEK8dOmSXL58+dBVrW6zzHqRRyMFq9Feizv6evrpp6Np+fPnz0cfhswHIXsxjytOd0RocjQ7OxuNRM1FXu7FPFaU9+\/fT5ymtyPQ7e3t0r\/bEavZh6TsZpFglmUQZaMlmdejIVBIUVpQ1c4x1iJKzcFgnZNLID5NacUXH2GamYikK1rdc96vvPKKvPTSS9GLtbMVcWm5U\/jxaVyXUvwq3Pg0LecoT26m2LP6IVBIUdYyoswy9Vo\/h589zUogLkr36mdbw8wcaEUZF5z93YxAv\/nNb0YjWHu+E1FmPWosBwE\/BAopymrnKC3a+IgyPtWa9PURP4eFKieNQNKFL\/FzkUaUmqnXpHPXSSNQd+q1Gp8sF+pkWYap12qk+fciECikKLMe2LgozXp5vh6SdTssd\/IJJEklfqWqPRed92Ie8+rNdK25eCz+Y89zpp1PjF8Zy9Tryc8Se1h\/BBBlhWOWJMq8Nxyov0iwx0kE0kZf7vlD96ItH18P6enpKbtQJy7LNEma\/c8yWsyyDCNK3g8QEKlbUSZ978wcUO6tSqwhAAEIQMAngboVZRKEalex+gRHLQhAAAIQKAaBhhFl0jRpMQ4hrxICEIAABI6SQEOI0k7Ddnd3V7z5+Pr6+lGypDYEIAABCGQg0N7enmGpk7NIQ4gyy2jSSNJ8yXtlZeXk0GdPIAABCBSQwBe+8AUxN96oF2E2hCiznJu0l+ybg3P69OkCRjP9JZsPDzMzM1FwYVPOCTbkRtMsyE313CwtLSFKTbg069hp176+vooPzbWirKeDo+GhWQc26dRgAxveUxoCjZWbuh9Rpj3BI36YaHiNFVy\/b13YaHjyniI3RclN3Ysy\/oistAPHmzo90g8ePJDXXntNXn75ZWlqatJkv2HXgQ250YSb3DTWh4i6F6W5u8jy8nLio4bcQ4Uo04P7ySefyO9\/\/3s5c+YMooxhgg250YiS3CBKTW6OfR1EScPThJCGR27IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNWWf\/Vb+dY\/X5F7\/zFRNw8LDYRGkAEy0GSN3JAbTW7qcdBS908PyXqgpu7+l0zdfSC\/+tcORBmDRsOj4WV9H7nLkRtyo8kNotRQC7QOouRNrYkaMiA35EZDgKlXv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLrcZq9js5pkxnZ6fMzc1JS0vLoaqIkoaniRqiJDfkRkMAUfqlVkO1tbU1GRoakmvXrklXV5csLi7K8vKyTE5OSnNzc1llREnD00QNUZIbcqMhgCj9UquhmhHjgwcPZHR0tGoVREnDqxqShAUQJbkhNxoCiNIvNWW13d1dGRsbk+7ubunt7a1axYryrW81R7ewa21trbpOURZABshAk3VyQ27y5GZzczNafH19Xfr7+2Vpaalubidat\/d6taJ87rnn5ObNm7K6uprpHOVjPxmODtbAwIAMDg7mOc4Nu+ze3p5sb29HHx6ampoa9nVqXhhs0qnBBjZ53lPz8\/OysLBQWgVR5qGnXNaK8uHDh6ULeKampmRjY6PiOcqbX9qR06dPR1JgVHkA37Dc2tqKPt0hyvJAwib9DQob2ORp32ZEaf6srKzIzMwMI8o88LTLJk29mot7RkZGZHp6Wjo6OspKc46SaSJN1pheJDfkRkMgfR2eHuKXZ9VqZgR59uzZ0jlKI8qrV6\/K9evXD31FBFHS8KoGKmEBREluyI2GAKL0S62GauaTiZGl\/e6k+X\/zk3QVLKKk4WmihijJDbnREECUfqnVWM294UBPT0\/i+UmzCURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VKrodru7q6MjY3JnTt3SlUuXrwoo6Ojh6oiShqeJmqIktyQGw0BROmXWg3VdnZ25NKlS3L58mXp6OioWAlR0vA0UUOU5IbcaAggSr\/Uaqi2trYmV69elevXr0tLSwuiVLJEBshAEx1yQ240ubl375709\/fL0tKStLe3a0oEX+fU\/v7+fvCtetqgAT41NSVzc3OIsgamNDwaniY+5IbcaHKDKDXUalhncXFRxsfHSxU6OztTpWmnXt\/6VnP0Kaa1tbWGLTfWqjQ8Gp4m0eSG3OTJzebmZrT4+vo6I8o84Gpd1owmNzY2ZHJyUpqbm6PRpfu7W9+K8rGfDEd\/PTAwIIODg7XuQkOsv7e3J9vb29GHh6ampoZ4Tb5eBGzSScIGNnneZ\/Pz87KwsFBahanXPPQ8LmvOWY6MjMj09PShi3usKG9+aUdOnz4dSYFR5QF8c\/Xw1tZWNNJGlOWBhE36GxQ2sMnTvs2I0vxZWVmRmZkZzlHmgedz2UoX93DVK9NEmqwxvUhuyI2GQPo6nKP0y7NiNfsdyu7ubunt7Y1GReY7lW1tbXyPMudxQAbIIGdkosXJDbnR5AZRaqjVsE78hgM9PT2l85XxsnZE+e7LX5QnWx6tYauNtyoNj4anSTW5ITea3CBKDbVA6yBK3tSaqCEDckNuNASYevVLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRLLVA1REnD00QNUZIbcqMhgCj9UgtUDVHS8DRRQ5TkhtxoCCBKv9QCVUOUNDxN1BAluSE3GgKI0i+1QNUQJQ1PEzVESW7IjYYAovRCbXd3V8bGxqSvr0+6urpqrrm2tiYjIyMyPT0tHR0dh+ohShqeJmSIktyQGw0BROmF2s7OjgwPD8vo6GjNorTSvX\/\/vty6dQtR5jxCyAAZ5IxMtDi5ITea3Ny7d0\/6+\/tlaWlJ2tvbNSWCr3Nqf39\/P\/hW\/7zBxcVFWV5elsnJSWlublbvhgE\/NTUVrc+IMj9GGh4NL39qEGUlZrynGFFq3lOH1rEjytXV1cR6nZ2dMjc3Jy0tLRW3Z+pcuXIlmsI1skSU+Q8Pb2pEmT81iBJRalIjwohSx62mtcyo1Pw8\/fTTmc5R\/vTbZ+TJlkeltbW1pu020sqIElFq8kxuyE2e3GxubkaLr6+vM\/WaB1yty5oLeObn5+V73\/teBD\/LxTyf\/tmoPPLxBzIwMCCDg4O17kJDrL+3tyfb29vRh4empqaGeE2+XgRs0knCBjZ53memVy8sLJRW4RxlDnpmRDg+Ph6tcfv2bXn\/\/fczn7c0U63PPPNMdDFQ1qteZ589JWc+czCiZFR5cKDMxVBbW1vRiXVEWR5e2KS\/mWEDmxytXsyI0vxZWVmRmZkZLubJCs+IbmNjIxoJvvjii9EVsObcpPnaSFtbW\/R72k+lc5xGuPGvnPD1EKaJsubSXY7pRXJDbjQE0tfhHGUOnu7XQ86dO1f2VRF7FWuWi3nsJrOOKN99+YvROUp+\/kIAGSADzfuB3JAbTW4QZQ5qiDIHrCNelIZHw9NEjNyQG01uEGVOavZ7lO7Uqx1dXrhwQXp7e3NWTF+cqVfe1JowIQNyQ240BJh69UrNfrpwi05MTHiVpKmNKGl4muAiSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTphZp71WrSVapeNuIUQZQ0PE2mECW5ITcaAojSGzV7M\/M7d+6U1cx6+7o8O4IoaXh58mKXRZTkhtxoCCBKv9Ri1aw8Hz58mOler1l3BlHS8LJmxV0OUZIbcqMhgCj9UvtzNfeiHi7mORLEqUWRATLQJI7ckBtNbvh6SE5q7tTrUUy3urvDiJI3dc54RosjA3JDbjQEGFF6oebzwc1ZdghR0vCy5CS+DKIkN+RGQwBReqVm7vc6OzsrPT09NT\/AudKOIUoania4iJLckBsNAUTpl9qfq5l7tQ4NDUU3ST+Kr4sgShqeJriIktyQGw0BROmXWkI1M8o0J3zz3BS92k4hShpetYwk\/TuiJDfkRkMAUXqj5o4i40V9X\/mKKGl4muAiSnJDbjQEEKUXau6deXxLMWkHESUNTxNcREluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAovVIzj9oaHx8vq3kUI0xEScPTBBdRkhtyoyGAKL1RM5J84403yi7asVOyPI\/SG+aqhZABMqgakoQFyA250eSGO\/PkoFbphgMGpLnylatecwCtYVEaHg1PEx9yQ240uUGUOaghyhywjnhRGh4NTxMxckNuNLlBlDmpMfWaE9gRLU7Do+FpokVuyI0mN4hSQY2LeRTQPK9Cw6PhaSJFbsiNJjeIUkOthnXiD36+ePGijI6OJlbkqlfe1JqoIQNyQ240BNLXQZR+eVatZi74MT9GjtWulkWUNLyqgUpYAFGSG3KjIYAoa6ZmnxhiCvl8DqUrzvhOIkoania4iJLckBsNAURZEzUjM\/OEkMnJSWlubpakC3o0G6j2fEtEScPT5ApRkhtyoyGAKNXU7DnFvr4+6erqiuok\/V3eDWR5pqUV5U+\/fUaebHlUWltb826mYZdHBshAE25yQ27y5GZzczNafH19Xfr7+2VpaUna29vzlDi2ZU\/t7+\/vh9p60qjPirK7u1t6e3tr2pX4aNUtZkX56Z+NyiMffyADAwMyODhY0\/YaZeW9vT3Z3t6OPjw0NTU1ysvy8jpgk44RNrDJ8yabn5+XhYWF0iqIMoXeUYvSPLZrZGREpqenpaOjo2wvrChnnz0lZz5zMKJkVHmAyHxY2draij7dIcry8MImvRXCBjZ5RGlGlObPysqKzMzMMKJMg3fUoqx06zvOUTJNlOdNbZdlepHckBsNgfR1+HpIFZ6+Rele5WqncNva2hK\/S4koaXiatzuiJDfkRkMAUaqpuQ9rrlYky9dG4jcc6OnpKV1NG6+PKGl41TKX9O+IktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pRaoGqKk4WmihijJDbnREECUfqkFqoYoaXiaqCFKckNuNAQQpV9qgaohShqeJmqIktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pRaoGqKk4WmihijJDbnREECUfqkFqoYoaXiaqCFKckNuNAQQpV9qgaohShqeJmqIktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pRaoGqKk4WmihijJDbnREECUfqkFqoYoaXiaqCFKckNuNAQQpV9qgaohShqeJmqIktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pRaoGqKk4WmihijJDbnREECUfqkFqoYoaXiaqCFKckNuNAQQpV9qgaohShqeJmqIktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pRaoGqKk4WmihijJDbnREECUfqkFqoYoaXiaqCFKckNuNAQQpV9qgaohShqeJmqIktyQGw0BROmXWqBqiJKGp4kaoiQ35EZDAFH6pVZDtZ2dHRkeHpbV1dWoSk9Pj0xOTkpzc\/OhqoiShqeJGqIkN+RGQwBR+qWmrLa7uytjY2PS3d0tvb29Yn9va2uT0dFRRJmDKzJABjniUlqU3JAbTW7u3bsn\/f39srS0JO3t7ZoSwdc5tb+\/vx98q0e0wcXFRVleXk4cVTKi5E2tiR0yIDfkRkOAEaVfah6rZRHlT799Rp5seVRaW1s9brm+SyEDZKBJMLkhN3lys7m5GS2+vr7OiDIPOJ\/L2vOVFy5ciKZi4z92RPnpn43KIx9\/IAMDAzI4OOhzF+q21t7enmxvb0cfHpqamur2dRzFjsMmnSpsYJPnPTc\/Py8LCwulVZh6zUPPw7L2\/KQpVe1intlnT8mZzxyMKBlVHsA3\/La2tqLzBYiyPJCwSX+DwgY2edq3GVGaPysrKzIzM8M5yjzwal02iyTNNjhHyTSRJmtML5IbcqMhkL4OF\/P45Vm1WrUrXd0CiJKGVzVQCQsgSnJDbjQEEKVfajVUm5qako2NjdTpVkSZDS4yQAbZklK+FLkhN5rcMKLUUFOuE7\/ZgC3T2dkpc3Nz0tLSUlaZESVvak3UkAG5ITcaAowo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVECUNTxM1REluyI2GAKL0Sy1QNURJw9NEDVGSG3KjIYAo\/VILVA1R0vA0UUOU5IbcaAggSr\/UAlVDlDQ8TdQQJbkhNxoCiNIvtUDVrCjf\/M7n5cufeyzQVutjM8gAGWiSSm7IjSY3fI9SQy3QOoiSN7UmasiA3JAbDQFGlH6pBaqGKGl4mqghSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTpl1qgaoiShqeJGqIkN+RGQwBR+qUWqBqipOFpooYoyQ250RBAlH6pBaqGKGl4mqghSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTpl1qgaoiShqeJGqIkN+RGQwBR+qUWqBqipOFpooYoyQ250RBAlH6pBaqGKGl4mqghSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTpl1qgaoiShqeJGqIkN+RGQwBR+qUWqBqipOFpooYoyQ250RBAlH6pBaqGKGl4mqghSnJDbjQEEKVfap6qTU1NydmzZ6W3tzexIqKk4WmihijJDbnREECUfql5qGYkOTs7KxMTE4hSwRMZIANFbITckBtNbngepYZaDevs7OzI8PCwPPbYwYOYv\/GNbyBKBU8aHg1PERtEWQEa7ylGlJr31JGsY0T50UcfSVtbm4yNjUl3d3dVUf77Pz0hX\/q7x6S1tfVI9qkei\/KmRpSa3JIbcpMnN5ubm9Hi6+vr0t\/fL0tLS9Le3p6nxLEte2p\/f3\/\/2LbuacO7u7uZRfk3b09L0we\/kYGBARkcHPS0B\/VdZm9vT7a3t6MPD01NTfX9YjzvPWzSgcIGNnnebvPz87KwsFBaBVHmoedh2TyiNCPKp\/7qT5EUGFUewDf8tra2ok93iLI8kLBJf4PCBjZ52rcZUZo\/KysrMjMzw4gyDzwfy+YR5Zvf+bx8+XMH5zT5OSDAFFp6EmADG02fIDfp1LiYR5MoD+sgytog8qZGBpoEkRtyo8kNotRQ87AOoqwNIg2PhqdJELkhN5rcIEoNNQ\/rIMraINLwaHiaBJEbcqPJDaLUUAu0Dnfm4U2tiRoyIDfkRkOAc5R+qQWqhihpeJqoIUpyQ240BBClX2qBqiFKGp4maoiS3JAbDQFE6ZdaoGqIkoaniRqiJDfkRkMAUfqlFqgaoqThaaKGKMkNudEQQJR+qQWqhihpeJqoIUpyQ240BBClX2qBqiFKGp4maoiS3JAbDQFE6ZdaoGqIkoaniRqiJDfkRkMAUfqlFqgaoqThaaKGKMkNudEQQJR+qQWqZkX5b33\/IH3\/yLMoXezIABlo3obkhtxocsOdeTTUAq2DKHlTa6KGDMgNudEQYETpl1qgaoiShqeJGqIkN+RGQwBR+qUWqBqipOFpooYoyQ250RBAlH6pBaqGKGl4mqghSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTpl1qgaoiShqeJGqIkN+RGQwBR+qUWqBqipOFpooYoyQ250RBAlH6pBaqGKGl4mqghSnJDbjQEEKVfaoGqIUoaniZqiJLckBsNAUTpl1qgalaUo8+dldHn\/jbQVutjM8gAGWiSSm7IjSY33JlHQy3QOoiSN7UmasiA3JAbDQFGlH6p1VhtcXFRxsfHoyoTExPS29ubWBFR0vA0UUOU5IbcaAggSr\/Uaqi2trYmIyMjMj09HVWx\/9\/R0XGoqhWluSG6uTE6P38h8ODBA3nttdfk+eefl\/b2dtA4BGCTHgfYwEbTLJh61VCrYR0zmlxeXpbJyUlpbm6WqakpOXv2bOKoElE21ie8GmKTa9V6fFPneoE1LAwb3lOa+NRjbk7t7+\/va17sSVjHiNH8jI6ORv+N\/+7uoxWl+bs3v\/N5ebLl0dSX8MjHH5yElxdsH9bX16W\/v19u377NiDJGHTbpMYQNbDRNyuZmaWmpbvpN3YvSHUGaEaaZDrLidA\/i27\/9H\/mXV38i\/\/nIU9FfN33wG\/nfz\/59aREjx\/\/7689qjjvrQAACEIBADgJnH\/2TvHO1J8cax7toYURpMJtPMg93PikRTxPjf3\/0l2WO9\/CwdQhAAAKNScBcL1IvP3UvSgM6y9RrvRwQ9hMCEIAABE4WgboWZXyqtdLFPCcLO3sDAQhAAAL1QqCuRZnn6yH1ckDYTwhAAAIQOFkE6lqUBmXWGw6cLOzsDQQgAAEI1AuBuhdlvYBmPyEAAQhAoD4JIMr6PG7sNQQgAAEIBCKAKAOBZjMQgAAEIFCfBAohSnM17OzsbHSEzN1nurq66vNo5djr3d1dGRsbk76+vrLXW+mc7s7OjgwPD8vq6qp0dnbK3NyctLS0lLZa7+eDzcVfQ0NDsrGxEb2mixcvlt2coshsDA\/3+Jvf4w8ZKDof+0Ywt2AzPcV9fxSZTTw38ew0ApuGF6Ub6vfee+9QwHO4p24WtZK8c+dO2QeDalcJu7cAjN8OsNq6Jx2OfTOb79yaD0r29wsXLkT3Bq72+hqZjTl2NjPd3d0lHuZDxbVr1yJeRedj821zY363oiw6G\/P6r169KtevXy\/7YG0YNQqbhhel2+DSRlknvcnn2T87ajp\/\/rw8fPgwGjHZEXSlm8jHRRIPf54b0OfZ3+Nc1s0GbMqPRFyc8DngYzi89dZb8oc\/\/KEkyqKzMYOR119\/vfRwCjdJjcKmoUUZf7PHfz\/OJn1U2\/7d734XlTZPUzHTqK4oK91E3v3kZx5TFv89zw3oj+q1+a5baZTo\/lsR2eQ5\/kXhY17n\/Py8fO1rX5MbN26URFn091Xp3thQAAAF9UlEQVRchu77tFHYFEKU7nm6oty9Jz5CNOGNv3b3zkbxEaRZ\/9KlS3L58mUx4qy0rm+BhahnR952ahE2B9TdaXv3HC58Dt4\/zzzzTMTJPUdZdDbuNSCGTU9PT+qjD+u15yDKEF35GLaBKNOhWzZmStq9T3Dak2iK9iHCFWZbW1vEqOgyMNOLv\/jFLyIW8Yt5iszGfrCyOYn\/3ihsCiFKe4FCEaZe4xcdMPVaLswkSdoRgvlv0g32izK1GP9o4Qrh5s2b0T8XkY\/pGz\/84Q9lcHAwml1JEmVR2SR9HG3E3DS0KOPTjUW4mKeSKCvdRD4+1Zp0MY\/7rM96nMKOX+kav+gg7fUVgU1Sw3PPPb355ptlz3p1j3+j84l\/rciyMqOoW7duyTvvvFNYNmmitBf3NEpuGl6URfx6iAlv0tRro1yqrZnJjk8JxWsUmY071WpnX\/j6THrK4iPKImenKF+7anhR2lFl0W44kCRKw6IRvvyrEWXaqMC98KCobOKzEOaGE+aHGzIkJ40bDiSfzmjk3BRClJrGyjoQgAAEIAABQwBRkgMIQAACEIBABQKIknhAAAIQgAAEECUZgAAEIAABCOgIMKLUcWMtCEAAAhAoCAFEWZADzcuU6Ivi\/f39iSiSHitWZGbudyMff\/zxQ\/cNtmzy3MSj0lMmisya137yCSDKk3+M2ENPBJIu6\/dUuuHKxG8oEL\/BvkaUZp34TS8aDhwvqCEJIMqGPKy8qCQCiDJbLpLubetLlPG7+GTbI5aCwPESQJTHy5+tBySQRZRmxPPrX\/9a\/vjHP4p58PXExET0IOP4tK39e7v77r+badwXXnhB7t69Gz1FYX19XUZGRmR6ejq6V6j5id8\/1n1qh\/l3dyrY3jzC1DT3XLVf7Hb3Ib6+vb2aqRXfdtrNKOxriT8aqdLy7tTrs88+G03R2v2z9ey+2NfOqDJg6NmUFwKI0gtGitQDgayiHB8fl9u3b5c98PqNN94oPX8wfns3U\/e73\/1udN9Pe9Nscy7U3vWnmijtOcALFy5EUrZTlHab5ncjIPMzNzcXPUXe3WZ7e7uMjY1F\/27EbJ5FamRk1n\/11VcjQdtb05llKnFIkmJWUdp9t1lIu21gluNQD3liH4tDAFEW51gX\/pWmXczjjnisYKyQ0m6kb5t9kogMaDMq29jYyDSi\/PDDDw89IT5ppOaK1JWXEW181Oge7KSnzJt\/t08CcZeNj3TNv9ltxUeK7nrxEXacgZG3\/UnaRuHDCYATTQBRnujDw875JJBlJBOXSiVJmOnRV155RV566aVIOub5lvbHbMs+QaHaiNI8fcKMYpN+jIDslKa7DVeUVkpW7vE6rpgqXcFq1qskyvhrNMunXfUa\/8Dh7lO1qV+fx5xaEPBBAFH6oEiNuiCgEaW9mfq1a9fKRGhfcKURZ1ZR\/vjHPy6NPt2Rl91GtenQaqK0\/24eTP3UU09Fo90sUrXnFPNOvcanouPhQJR18XZhJx0CiJI4FIaARpSVnmFpwcUvfjF\/745Mk0aUrkzMiNI9B5pFLHmmXk09O8L91Kc+JeZP0rSrjxFltQ8WadsoTAh5oXVJAFHW5WFjpzUENKK00rtx40bpYh07QnPPQQ4NDYkddVpZnD9\/PjpHaUad5mIcMzVrBGVHoffv349q2ulQ++9WJrbmuXPnDn3h3xWlmQI2F\/OYc61WgPEpVHcK2b1QKc4xaYScdURpp4jd15F0nLIcB83xZR0IHBUBRHlUZKl74ghkadDxc5T2RbjPqjR\/5z7H0hWbkacR1le\/+tXoKyb2KlT3eZjm3y9fvhx91cN+ZSTpXKgVWrWpVyOmtK+H2OlTK3fDIG3a1X2tDx48KEk3qyjN+mnnWl05J43AT1xY2CEIOAQQJXGAwBEQSBPuEWwqc8msgjrKW81xw4HMh4sFTxABRHmCDga70jgETpoo8wrKvYWdz6PCzQZ80qRWKAKIMhRptlMoAidJlHbaOOm7jmkHJa9YsxzcoxypZtk+y0BASwBRasmxHgQgAAEIFILA\/wOo5rXU7EfjPwAAAABJRU5ErkJggg==","height":275,"width":458}}
%---
%[output:861b919d]
%   data: {"dataType":"textualVariable","outputData":{"name":"tau","value":"0.0159"}}
%---
%[output:5ca1803f]
%   data: {"dataType":"textualVariable","outputData":{"name":"num","value":"1"}}
%---
%[output:1ad2475a]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"den","rows":1,"type":"double","value":[["0.0159","1.0000"]]}}
%---
