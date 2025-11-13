clc ; clear;
%%

%% ========================================================================
% 1. SYSTEM PARAMETERS
% =========================================================================
fprintf('Defining system parameters...\n'); %[output:84378190]
M_w = 0.5;   % Mass of one wheel
M_r = 1.0;   % Mass of robot body
m_b = 0.02;   % Mass of the ball
R = 0.05;    % Radius of the wheel
h = 0.1;     % Height of robot CoM from ball platform
r = 0.002;    % Radius of the ball
g = 9.81;    % Gravity
d = 0.03; 

% Moments of Inertia
I_w = 0.5 * M_w * R^2;        % Inertia of one wheel
I_r_com = (1/12) * M_r * (h^2 + d^2);       % body about its own COM
I_r = I_r_com + M_r * d^2;              % Inertia of robot body (given value)
I_b = 0.4 * m_b * r^2;      % Inertia of a solid sphere ball

% Derived parameters
l = h + r; % Height of ball CoM from wheel center
%%

% --- 1. Define Parameters ---
M_w = 0.5;   % Mass of one wheel
M_r = 1.0;   % Mass of robot body
m_b = 0.02;  % Mass of the ball
R = 0.05;    % Radius of the wheel
h = 0.1;     % Height of robot CoM from ball platform
r = 0.002;   % Radius of the ball
g = 9.81;    % Gravity
d = 0.03; 

% --- 2. Moments of Inertia ---
I_w = 0.5 * M_w * R^2;        % Inertia of one wheel (Note: Not used in M1, M2, M3)
I_r_com = (1/12) * M_r * (h^2 + d^2);       % body about its own COM
I_r = I_r_com + M_r * d^2;              % Inertia of robot body (parallel axis theorem)
I_b = 0.4 * m_b * r^2;      % Inertia of a solid sphere ball

% --- 3. Derived Parameters ---
l = h + r; % Height of ball CoM from wheel center

% --- 4. Define Matrices ---
% M1 = [[(I_b/(r^2))+m_b , -m_b*l],
%       [-m_b*l,            I_r+m_b*(l^2)]]
M1 = [(I_b/(r^2)) + m_b, -m_b*l;
      -m_b*l,            I_r + m_b*(l^2)];

% M2 = [[1,  1+(g*m_b)],
%       [1-(m_b*g), 1 - (M_r*d*g)-(m_b*g*l)]]
M2 = [1,             1 + (g*m_b);
      1 - (m_b*g), 1 - (M_r*d*g) - (m_b*g*l)];

% M3 = [[m_b*l],
%       [-m_b]]
M3 = [m_b*l;
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
C = eye(4) %[output:717ff987]
D = zeros(4,1) %[output:052d9a1f]
%%
[num,den] = ss2tf(A,B,C,D) %[output:17d34675] %[output:9a8921cb]
%%
% --- 2. Create the individual Transfer Functions (TFs) ---

% H1(s) = Y1/U (from num row 1)
H1 = tf(num(1,:), den) %[output:8546a222]

% H2(s) = Y2/U (from num row 2)
H2 = tf(num(2,:), den) %[output:2e33cd79]

% H3(s) = Y3/U (from num row 3)
H3 = tf(num(3,:), den) %[output:0f51ec1d]

% H4(s) = Y4/U (from num row 4)
H4 = tf(num(4,:), den) %[output:88dfa763]
%%
sys = ss(A, B, C, D);
%%
rank(ctrb(A,B)) %[output:07841c68]
rank(obsv(A,C)) %[output:282facca]
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
Q = [100 0 0 0 ; 
    0 10 0 0;
    0 0 1000 0;
    0 0 0 10];
%%
K = lqr(A, B, Q, R);
%%
% Closed-loop system
Acl = A - B*K;
sys_cl = ss(Acl, [], eye(4), []);
x0 = [0.005, 0, deg2rad(88),0] %[output:13ba49e1]
x_ref = [0.000, 0, deg2rad(90),0] %[output:256dedc6]
x0 = x0 - x_ref;
% Simulate response
t_sim = 0:0.01:20;
[y, t] = initial(sys_cl, x0, t_sim);

% Plot
figure; %[output:67ec4bae]
plot(t, y); %[output:67ec4bae]
xlabel('Time (s)'); %[output:67ec4bae]
ylabel('States [x, xdot, theta, thetadot]'); %[output:67ec4bae]
legend('x (ball pos)', 'xdot', 'theta', 'thetadot'); %[output:67ec4bae]
title('Closed-Loop Response of Ball balancer robot with LQR'); %[output:67ec4bae]
grid on; %[output:67ec4bae]

% Eigenvalues
disp('Closed-loop eigenvalues:'); %[output:9f4e2fbc]
disp(eig(Acl)); %[output:2cb5c6d2]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:84378190]
%   data: {"dataType":"text","outputData":{"text":"Defining system parameters...\n","truncated":false}}
%---
%[output:60988226]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["69.9103","0","72.8678","0"],["0","0","0","1.0000"],["469.3566","0","413.7732","0"]]}}
%---
%[output:06e878d2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["0.7015"],["0"],["10.6283"]]}}
%---
%[output:717ff987]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":4,"type":"double","value":[["1","0","0","0"],["0","1","0","0"],["0","0","1","0"],["0","0","0","1"]]}}
%---
%[output:052d9a1f]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"D","rows":4,"type":"double","value":[["0"],["0"],["0"],["0"]]}}
%---
%[output:17d34675]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":4,"type":"double","value":[["0","0","0.7015","0","484.2024"],["0","0.7015","0","484.2024","0"],["0","0","10.6283","0","-413.7780"],["0","10.6283","-0.0000","-413.7780","0"]]}}
%---
%[output:9a8921cb]
%   data: {"dataType":"matrix","outputData":{"columns":5,"exponent":"3","name":"den","rows":1,"type":"double","value":[["0.0010","0.0000","-0.4837","0.0000","-5.2740"]]}}
%---
%[output:8546a222]
%   data: {"dataType":"text","outputData":{"text":"\nH1 =\n \n                  0.7015 s^2 + 484.2\n  ---------------------------------------------------\n  s^4 + 3.22e-15 s^3 - 483.7 s^2 + 2.026e-13 s - 5274\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 48 46 55 48 49 53 32 48 32 52 56 52 46 50 48 50 52 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 51 46 50 49 57 54 101 45 49 53 32 45 52 56 51 46 54 56 51 53 32 50 46 48 50 53 53 101 45 49 51 32 45 53 46 50 55 52 48 101 43 48 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:2e33cd79]
%   data: {"dataType":"text","outputData":{"text":"\nH2 =\n \n                 0.7015 s^3 + 484.2 s\n  ---------------------------------------------------\n  s^4 + 3.22e-15 s^3 - 483.7 s^2 + 2.026e-13 s - 5274\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 46 55 48 49 53 32 48 32 52 56 52 46 50 48 50 52 32 48 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 51 46 50 49 57 54 101 45 49 53 32 45 52 56 51 46 54 56 51 53 32 50 46 48 50 53 53 101 45 49 51 32 45 53 46 50 55 52 48 101 43 48 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:0f51ec1d]
%   data: {"dataType":"text","outputData":{"text":"\nH3 =\n \n                   10.63 s^2 - 413.8\n  ---------------------------------------------------\n  s^4 + 3.22e-15 s^3 - 483.7 s^2 + 2.026e-13 s - 5274\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 49 48 46 54 50 56 51 32 48 32 45 52 49 51 46 55 55 56 48 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 51 46 50 49 57 54 101 45 49 53 32 45 52 56 51 46 54 56 51 53 32 50 46 48 50 53 53 101 45 49 51 32 45 53 46 50 55 52 48 101 43 48 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:88dfa763]
%   data: {"dataType":"text","outputData":{"text":"\nH4 =\n \n          10.63 s^3 - 9.44e-15 s^2 - 413.8 s\n  ---------------------------------------------------\n  s^4 + 3.22e-15 s^3 - 483.7 s^2 + 2.026e-13 s - 5274\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 49 48 46 54 50 56 51 32 45 57 46 52 51 57 56 101 45 49 53 32 45 52 49 51 46 55 55 56 48 32 48 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 51 46 50 49 57 54 101 45 49 53 32 45 52 56 51 46 54 56 51 53 32 50 46 48 50 53 53 101 45 49 51 32 45 53 46 50 55 52 48 101 43 48 51 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:07841c68]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:282facca]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:13ba49e1]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"x0","rows":1,"type":"double","value":[["0.0050","0","1.5359","0"]]}}
%---
%[output:256dedc6]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"x_ref","rows":1,"type":"double","value":[["0","0","1.5708","0"]]}}
%---
%[output:67ec4bae]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEeCAYAAAAgg6RKAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQ2MXtWZ34\/ThDKL2JhxDVmvjcekduhGi1mHVtNJWbNxrbDd2Juk3Z2PIrHTCZ0GvEmDzXjGJHZAiWcYhmRloNEIJlNLdGxKJTc2Vau63kCJRq6EkSbqthsjhjGMvPUSPiQ2dZSSTPVcc17OHN+P83Wf9973\/F8JMeM55zz3\/p7nPv\/3nHs+ViwtLS0JfEAABEAABEAABEolsAKCWypfNA4CIAACIAACCQEILgIBBEAABEAABBgIQHAZIMMECIAACIAACEBwEQMgAAIgAAIgwECgpQX35ZdfFv39\/eL8+fMNlKOjo6K7u7vx+0MPPSQmJyfFzMyM6OzsZECebuL06dOir69PDA4Oir1796YWqsq1ZkGS15f29x07doixsTHR1tbWNMZ1M6zyTIuLULzfeustMTAwkOCZmppK\/q\/+3t7efhm6ELEYoo26+ZSut6r3\/fTTT4uRkRGh50iVsUkZ3Sc2ddJyth77aWVUm83O5Xkx2bKCK51clPyrEvytLrjkh7wvE3VMnGVes4wHacNWcG14Q3DL9OTlbYfMORQnzz\/\/fOaXdJs7SxNGErfDhw+L+++\/P\/mybCOe0rZpnbycvWbNGjE9PS02btwoigRXLWtz\/xxlW1JwVYeo33ZUh8p\/Dxn8Pg5rJcHVv2FKf6xevTrpQaX1mHzYtWJdn3iw5Q3B5Y2gUDlH5rOyvsjKGFRHp0zF07ZXrOZstYd98eJFMTw8LE6cOCHkdSwuLiYjl5\/61KcuGzWTbPN66bzeXm6tJQU3LyjIIR0dHY1h5bTg179Bbd68eZlQqEGQ1QPRv63pIqT2YKj9u+66S+zatSvIkHLR9dM1F5WRXB577DHxxBNPiLm5ueRW84ZrshKJTOhvvPFG41uq\/DfZrj7knPYtVtpWBUJyo2vT\/UT\/pvtBTU5qOzSMv2fPnsbrB\/U+866FbOjxkHYd+kOex1+\/5qxv7Fm8ZdtqQtJ7zKovfQRXj4+sZCnvX\/Vz2vXrw+QqS\/W5XlhYSF4F0acodvS\/58Weyu7qq69O4j3tdYi8lnvvvVecOnUqeT7kvbs+W7pI+MaIzle9bxnf6v3ecsst4sCBA8l9rF+\/PnnFpX6ozrlz5xrDzvQ3GoIuGk0xEem8Lw96\/iB7WYJb9pcQX7FuScG1+Qapl01LTARZJr21a9c2vnHp8OUDkzU0IoM8y0ZR4JrcV9H105CMSZm894NZopt1fdKeTJ50n\/SOUIqtnozlN1j13bvqg1WrVqXW10U36x5kAiWhTLsOtZ0333zzsnkAJvGQJ7pF\/F966aVGIlNtke\/UT56PTIbgdH9Q27bvcNMSkD56lPWc6PGS9dzIL0l5Q45ZYqfHVpbP9d6TGntpPci0a6H7po8uVLoPTZ6tEDGij5KobUpeahnqiMh3uEWCe+utt4oXXnhhmWuzepVFgqt+YTXJLfT8pwmu2g56uL5fDQzrmzgvLWmRoykByeEL9SGTDwj921e+8pWkzJkzZxq9NbW9tN6c3uM4dOhQ8u1c2lCv2WfSVFY7addPQzRZ90i9PVknrUeSNQEqL5EQI\/0LSZp98kPWAyU5q9\/UZZv6v23ZsqUhlHrPmISe\/m3Tpk0NwdUTNtmid0Z536bpb7bDbiY+Iv42Q8pZj0aW6OsxKr\/AuAhu2nBjVnzovY+iL5A62+PHjyeCoN6XGttq3OrPlnxe5ZeZotgjwc0bzZH3ol6LqW\/Tni3ZHrHbv3+\/ePDBB5Nh1LxntKg3l5V3yM\/SRzIX6b1XmliaF9tpIw9ZuSuE4KptqM92WuybjDAZyknwYi3Xw\/URXJmAZeKR7xr1d2InT57M7IHkvdCnQKAhuPHx8csEW02wUtTpgZMfeiBockTejGp9aDDt+h9++GFx3333Jc2q71P1e6RhZN1W0bvBLMHVh0TzhJmEb+fOnZeNIqhJPOs+1QS0devW1FnfahkajqYerjrULeNHJui0EY00kUl7MtMSkImPyC9nz551nrWeNnRI16cPpUq\/+Ahu2tC7\/q4+a1i\/aERE751KwVW5pj03WV+Gqb2i2JPJvGi+QZrYmfq26NkyfUZlHsoSOjWWH3\/8cfH9739fvPbaa+L6669P\/j8xMZG8RpG5QLan937T4j3LB2krLEIIrvpuNk9wy3qfHUp5W05wCUzeN7+8d7imgiuFLC2RfPGLX0wdgqTrikVw04bO1SGeoqQnl22lvStXe6b6lwZXwVXb0QVXDuNmXYt8p1UlwVWFhbjLBCXf6cphVflFoyzB1YWlqIerJ+asHq4aS2UJbtqEHNXHaSJSNcFVcyENcz\/33HNJz1YOHZM4qSNtRfyzZikXjcYUCW5azpbPW1dXl9i+fXtjJCptBIzqy5HJKvdu6TpbUnCzZimr7zDSZimbDCmnfYNTe37y26naa9KTsT4MVpchZXXYK21NbVqPJY150VBYmnil9UzViSr6+xubIWUTwU1LtpSwZE\/adJ2x6bBjURJTRVUf+tR7uHJkJGvY3Edws4aUv\/SlL4kvf\/nLQu0pZg3\/ql+i1OdGjzfZw80S3LwhZRotUodNs3pCaRPO8uIxa5KY7esalyFlkxjRR9xUwZLvqWX86MKYxiJNPIuuw0Rw02Yp6++x5YgM+UN\/h5s2m7mKa\/5bUnDVb0xpD0veTMmiyQp5E3aKJnfIv\/tOmkq7J\/nNTg5F6mXUYd2ie6RencnEjqwvEroAyLaKJk3Ja6SJSmkTT8ie\/u5Vv4a093t6GX3SVJ7gFl2L+iVNt5P1DtCEf1ESKxoelSMqNDytvwKR1xliSDktFot8JJ8D9QtaHkfpLxPBzXqlYxp7acncVHCpnIlvTZ4tk3bUMlkz2VUhSpv4mTULnEaZ0sRanaUsR6KKYjVvspt63XnliK2Mm6wvRVlLi9L816x\/a1nBJaBpD5\/JTlNF0\/rThhf1Hk7WeyvpaPVhsV0WlCe4NNxddP1pbPShGJkUnnzySfGDH\/wgmcAhBS9rR66sd3Jqj0s+NFnvE+UQbp7vspbz+CwLku+z04aUi+Ioa7g5b+eyIh8VJbEiwc2azEP16Bmgj5yRKofs5BcP+r\/pTlN5y8b0GNffGepDzioTSsTyvaM+4Smvh2sS23mx59PDlc9lkW\/lc1K05K6oHf0+sr7gpY1M6aMNdO16T1SPaz1uQgtuVs6WXE3W4cp7qOrmF5UVXFWwiqZ46w82Nlfw\/\/5WNIPU34J7C1nvytxbRE0QAIEqE1Df6apb81b5mtOurZKCS9\/shoaGktm89JE\/6+sQ1W9E9O2ZehQk1LOzs9i31zMSIbieAFEdBEAABDQClRRcXTT1mcXqPVBZ2nUma8N\/eNyNAATXjRtqgQAIgEAWgUoKLiV7+kgR1X+XN9MqwwwITxAAARAAgdYnUFnBVfc7zurFSsH97Gc\/29jv12Qd1g033ND6nsUdggAIgECLEqD9qzds2FC7u2sJwaVdU+REKeoN0\/qyvLNXSXDn5+dr56y6XjB483oOvHl5kzUw52VeV96VFVzXIWV1wlXaJCs8HLwPBniDNz8Bfot1FQB+UmEs1pV3JQVXH0LOmzSl\/40E9+DBg+KRRx7JPHe1rs4KE6r8rYA3L3Pw5uWNL5XgbUqgkoJrsyyI1uCS6KpDymrvOA0EEpJpeIQpB95hOJq2At6mpMKVA\/NwLE1aqivvSgouAc\/a+EJOlOrt7U3W3dJH3fjCZE\/bujrLJBCrWObVV1+t5QSHKrI0uSbwNqEUtgyYh+VZ1Fpdc3hlBbcIuM\/f6+osn3tuZl0kI1764M3Lm6yBOS\/zuuZwCC5vnERpDcmI1+3gzcsbgsvPG4LLz9zZYl2d5XzDTa4IAeB1AHjz8obg8vOuaw5HD5c\/VqKzCAHgdTl48\/KG4PLzhuDyM3e2WFdnOd9wkytCAHgdAN68vCG4\/LzrmsPRw+WPlegsQgB4XR4zb2zbyhtroa2Z7gAIwQ1NvsT26uqsEpGU2nTMAlAq2IzGY+aNZ7sZERfGpo3vbMqGubowraCHG4YjWskhELMANCMwYuZd10TcjDipmk0b39mUrdJ9QnCr5I0WvZaYBaAZLo2Zd10TcTPipGo2bXxnU7ZK9wnBrZI3WvRaYhaAZrg0Zt51TcQucUJb4B4+fFh87WtfEw8++KBQd98zaY926Dty5IjYv3+\/U30TG2oZsvf88883zjnX69v4zqas7XWWWd5ZcN966y0xMDAg5ubmrK6Pzqs9duyYVZ3QhevqrNAcuNqLWQC4GKt2YuYdy7NNW9x++9vfFnfeeadYu3atGB4errzgUozSvvdbt25tbMurxq2N72zKNuMZzLLpLbh79+5NhZdmUB40AMGtUgiUfy0xC0D5dC+3EDPvuiZi6UX19LO0feNlObW3KMvdfPPNySEudB746Oio6O7uToqre83T7zMzM0nOLurhFrVLPez+\/v7EHnWk5AEyst6JEycS++q1yF75\/fffL9ra2pYFr43vbMo24xmE4CoE6uqsKgWOzbXELAA2nEKVjZl33Z9ttef65ptvZg7Bqj1FKXCvvfZaInpUb2hoSIyPjychJX+m88HpUJjZ2VkxNjaWjE7mDSnntbtq1apkhLOnpycRdroeEl5q9\/jx42JhYSEZOqaR0AceeEAcOHAgOS5VvT\/9vHIb39mUDfVchWjHuYerGpfDy2m9Xf34vBAX7dtGXZ3le9\/Nqh+zADSDecy8W+HZpl7gt771LfGhD30o9VxvXbSkMHZ1dTV6tVlniMterY3gprW7ZcuWZeeOq0eqvvTSSw3BTYv\/rGFlG9\/ZlG3GMxi8hwvBrZIbq30tMQtAMzwTM289Ebff+8NmuMDY5lvf+b3UslmCSYWpg7N7926xb98+Qb3EtKFntT79PDk52bAjjzA17eGqk7Fku+vXr192Drl+TapNOYQtLwCCaxweywuqZ9bmNTE4OJg5M83RtFe1kY\/9N\/FPdv+2uPW+67zaQWUzAjELgBmhsKVi5l3Xno8aAdQLffbZZ5N\/oklR+tBrUQ9X7fHqwujTw1XbzevhqterCzHdEwTX83nPG1L2bLqU6hDcUrBmNhqzAPCSvmQtZt51F1z1nSe9i6VlP2kTjNLe4ZLvaah4cXGx8d6W2qCy9G6XJinRbGZZzrSHm9Zu3jvcQ4cOiY6OjmR4G+9wl2eAIO9wm5FUfGxCcH3o2deNWQDsafnXiJl33QX3q1\/96rLlPTSKKCcg6b1guaY1bTaxHMZVZwyvWbMmGYZ+5plnknfDZ8+eNZo0RfXkkLQ6PJw1S1lfMqrXyfoSYeM7m7L+T1S4FoIKbtoQszolPNxl+7UEwfXjZ1s7ZgGwZRWifMy865qIbf2eN9vXtq2s8nnLklxtYB3u0tKSKzy1Hont0aNHG2ux6G\/ym46cOh7CTog2ILghKJq3EbMAmFMKVzJm3rEILkVL3prWENEUWnCx05QQQXq4dVsWBMEN8TiatxGzAJhTClcyZt4xCW64iKlGSza+sylbjbu7dBUQ3Cp5o0WvJWYBaIZLY+Zd10TcjDipmk0b39mUrdJ9BhFcuqHQQ8rq++C898D6NmJ0LUXLkNDD5Q3BmAWAl\/QlazHzrmsibkacVM2mje9sylbpPoMJrhTdkZGRZffnMmlK3bWEGlO3JtPhpa3zKgIMwS0iFPbvMQtAWJJmrcXMu66J2MyzH5Qy3cEvbwMNW5tll7fxnU3Zsq\/bpv2ggmtjOK+sut8nrR3LCxoS54MHD6ZugZZlA4IbylNm7cQsAGaEwpaKmXddE7FtBEBwbxDz8\/O22JpevpKCSwJLH9qbmT767yo108BT65DgfvqPPyE+d2hd0x0QwwXELADN8G\/MvOsuuHmnBanrXum1GeU+eUKPeiKQ3LqRDhGQI44uI43csWvjO5uy3PeRZ89ZcG3Ow1WPbjK5eb1Hm7X4m9rS1\/6a2ILgmnghXJmYBSAcRfOWYuZd10QsvZt1WpC+xJJypBRc2k2KjsmbmJhIjsmj3aRoswrqsGBI2fy54SjpLLj6xeVNmrI5M1f2aOXWYFJU03ZbkWXlsVBy+Fn9PQ0iCe5f\/a3\/Kf773z6U\/PnUqVMcrKO1QVvN0SHZ+PAQiJn3tm3bajnUqEZG2mlB+qsz9Xe5YxRt60g5UB31e+KJJxrbLPJEn7sV+rKUlYvJr\/on2iHl0OtwbYaUdSeoE670Tb9lWfRw3R8Kl5ox97hcePnWiZm33sOd\/6MVvjhLrX\/DM+n7Duk9U\/3VmSq4J0+ebJxxS4Kr\/q1ugmsqonUdyQjSww0tuPoQss2wiMkkKghuqTnkssZjFgBe0pesxcy7rolYjZO004Ji6eFCcA0zRsh1uKbLgvSDl+Xv8v1F1qVDcA2dGqhYzAIQCKFVMzHzrrvgZp0WRAFA72blYfB4h4tZysm7g76+vmXJQT982DRzZG18oe\/vqW98IWfo0dAKBNeUdrnlYhaAcsmmtx4z77oLbt5pQeos5XvvvTcZOj5w4IBob29flnvVHCjzKGYpN+NJvNxmkCHlatyK+VWgh2vOKkTJmAUgBD\/bNmLmXXfBtfV1K5W38Z1N2SoxguBWyRstei0xC0AzXBoz77om4mbESdVs2vjOpmyV7jOY4KrDHfoNmqyN5YSCHi4n7bgn8fCSvmQNglu\/HYiaESdVs2kjojZlq3SfQQRXnby0c+fO5OV+b2+v2LRpkxgYGEgWYHd2dlbmviG4vK6IWQB4SUNw65qImxEnVbNp4zubslW6zyCCqy8LUpfx0ESqI0eOCLkouwo3D8Hl9QIEF7y5CNQ1EXPxqbIdG9\/ZlK3SPZciuOo6Wpe9jssGBMEtm\/Dy9iG44M1FoK6JWPJRTz+j3dnkaKHNCKHJXgRc\/rCxY+M7m7I211B22SCCSxep7g6liqy+C0rZN2TSPgTXhFK4MhDccCxNWoqZd10TMQRXCBvf2ZQ1eWa4ygQTXH0TChLgycnJZBPt6elpkbXNIteNqnZIcD\/5DzeIO459vBnmo7MZswA0w9kx865rIpZxQutwT5w4keTN733ve+LJJ58UN998c3IqEO0Rr66nVSeqyomp1A7Nm5mbm0sOMqB6tNeyuj+C694IZceyje9sypZ93TbtBxNcG6PNLgvB5fVAzALAS\/qStZh51zUR5\/Vw6W80B4ZEdM+ePUkHZtWqVcsmpFIHRx7aQodXyDPC6SShoaEhMT4+nnR69LPGmxGfWTZtfGdTtkr3GERwQ++lXDYgCG7ZhJe3H7MA8JKG4NY1EecJrtzOUX2\/S0JKIivPw1Xf29LfpODSLlTqp4qTWOX12fjOpmwznsEsmxDcKnmjRa8Fgsvr2Jh564n4reMbeOFbWmvf+eqyGnmTpnTB1bfRla\/vqEFVcOXrPWnIZPtby9sIUtxGRG3KBrm4QI14Ca5++HvWNQ0ODiZrcavyQQ+X1xMxCwAvafRw65qI83q4tKcBzVLWBTdruaV+Vq7aE0YPtxlP5Ac2vQRXDZIqbnCRhRaCyxt0EFzw5iIQi+Dq73DV09rUIWWaMCUFlw50oWVG9KnSvggYUuZ6OppkB4LLCx6CC95cBOouuHK1x5kzZxqzlNN6uDQBSp2lrK4GkXNqiPljjz2WTJiSM5\/37dsnnnnmGfHII48kpwxV6WPjO5uyVbrHID3cKt2QybVAcE0ohSsDwQ3H0qSlmHnXNRGb+LXVy9j4zqZslbgFFVz1nS6t9Tp37pyYnZ2t3PAFBJc3BGMWAF7SeIdb10TcjDipmk0b39mUrdJ9BhNcuQ6M1nzt2rUrmSRFC6\/pnQENd1Rp0tTB634sru+6ChtfMEUiBJcJ9PtmYuZd10TMGyHVtGbjO5uyVbrbIIKrrsPVTwiq4l7KEFzeEIxZAHhJo4db10TcjDipmk0b39mUrdJ9QnCr5I0WvRYILq9jY+Zd10TMGyHVtGbjO5uyVbrbIIJLNyS3DFOHlGVvt6enR3R3d1fmvtHD5XVFzALASxo93Lom4mbESdVs2vjOpmyV7jOY4NJN0fCxvvuJutl2VW6cBPej664Q97x4Y1UuqaWvA4LL696Yedc1EcsIaebxfCb7LKvXF\/pAGhvf2ZTlffryrQUV3CrdWN61QHB5PRWzAPCSRg+3rokYgovj+ZqRKxo21SVGpr1kWgiunoyRdQMQXF7XQnDBm4tA3QWX+3g+OSpJK0luu+028e677ybLOOlDK0xowwz60DJPuepEbqIR+thVG9\/ZlOWKPRM7wXq46q4numF5LqPpziaqcFJbJiKq7tBSFAgQXJPQCFcGghuOpUlLMfOuayLO6+HS38o4no+O8evv7xcTExMNMZW2Dh061DjuTz8WcPfu3YJ2rMKQssnTuLxMEMGVYhdqva3+LoHW+HZ0dOROvJLLj+j25NmP6OHaB0QZNWIWgDJ4FrUZM29dcOnLdZU\/+y7ctOzy0t7hlnU8HwmpegCCPNhg\/\/794sEHHxRyS0mZ3+k6tm\/fLiC47hEVRHDzzsN1uTQSWPrIzTL03\/U2yf4DDzyQBAiVheC6UC+vTswCUB7V7JZj5t2KPdys04J8j+c7fvz4sp0ApeCqK03olCL6yE4PBNfviQ4iuPIbkAwMv0v6wLlyKRH1eBcWFjJ3q6K\/02fLli1Gw88YUvb1kF39mAXAjlSY0jHzjklwfY\/nQw83zPNm00oQwSWDIXeU0oeQ8wSX3vcePnxY3H\/\/\/YLeSZi87yXBfXfFX4unf+1fJqxOnTplwwxlLQmQX9auXWtZC8VdCcTMe9u2bWJ+ft4VXdPrmR5AH+J4PuooyWNV5YQoAkDvi5v1DjcrF5Nf9U8d\/ewsuHIYmb4lFX1sJ03ZDClT2a1btyYHNGOWcpEnmvP3mHtczSAeM++693C5j+dT9074xje+IV555ZWk80IffZYy5Vibyam2sW\/jO5uyttdRZnlnwS3zovQebdakqTzRp2ns8v2Dfq0YUi7Te5e3HbMA8JK+ZC1m3nVNxM2Ik6rZtPGdTdkq3WcQwc2bNOUy1OyyLIigoodbpdD64FpiFoBmeCRm3nVNxM2Ik6rZtPGdTdkq3WclBZcAZW18kTdBy0ZwyYY+Jb9Kjmmla4lZAJrhx5h51zURNyNOqmbTxnc2Zat0n16Cq4pi3k0NDg5W7jxcCC5fGMYsAHyUMaJABOqaiJsRJ1WzaeM7m7JVuk8vwZU3EnodbtmA5GJ49HDLJn2pfQguD2dpJWbelIjxqS8B05nHUQtu3dwLweX1WMwCwEsaX3CawRtfKvmpQ3D5mTtbhOA6o3OqCMF1wuZcCbyd0TlXBHNndE4VIbhO2JpTCYLLyx3JCLx5CfBbQ4zzMofg8vL2sgbB9cJnXRnJyBqZVwXw9sLnVBnMnbA5V4LgOqPjrwjB5WWOZATevAT4rSHGeZlDcDN4q1uH0fF9RWfVcrgNgstB+QMbSEbgzUuA3xpinJc5BJeXt5c1CK4XPuvKSEbWyLwqgLcXPqfKYO6EzbkSBNcZHX9FCC4vcyQj8OYlwG8NMc7LHILLy9vLGgTXC591ZSQja2ReFcDbC59TZTB3wuZcKXrBpX2M+\/v7xfnz5y+DaHs8n7MXDCtCcA1BBSqGZBQIpGEz4G0IKmAxMA8I06CpqAVXHijQ1dUldu7cmZyj2NvbKzZt2tQ44DjrqDwDtsGLQHCDI81tEMkIvHkJ8FtDjPMyj1pw9b2U1fNraZbykSNHxNjYmGhra+P1SoY1CC6vG5CMwJuXAL81xDgvcwjuwEByIhD1ZNUD5F3Owy3bdRDcsgkvbx\/JCLx5CfBbQ4zzMo9acAk19WrpQ6KriuzJkyfF7Owseri88Vgpa0hGvO4Ab17eZA3MeZlHL7jqe9zu7u5EgCcnJ0VVNrtQwwE9XN6HA8kIvHkJ8FtDjPMyj15weXH7WYPg+vGzrY1kZEvMrzx4+\/FzqQ3mLtTc60QtuHkH0OMdrntQtUpNJCNeT4I3L28MKfPzhuAqk6ZU\/BBc\/mCsmkUIAK9HwJuXNwSXn3eUgkuzkUdGRgppDw4OJpOpqvLBkDKvJyAA4M1LgN8aYpyXeZSCKxHnDSnzusHMGgTXjFOoUkhGoUiatQPeZpxClgLzkDSL24pacIvxVKsEBJfXH0hG4M1LgN8aYpyXOQRXiGTDCznEPDMzI86dO+e8Bldta3R0VNBSo7SPXI504sSJ5M8mw9cQXN6HA8kIvHkJ8FtDjPMyj15wad0tHVwwNDQkdu3albyzpUMLaF9lWotr8w6XDkKgdsbHxxMvyp83btx4mVfVDTfk0HZPT0+mQFMDEFzehwPJCLx5CfBbQ4zzMo9acNV3uPqBBS6zlKl3q+5Ope7NXORWVYCzykJwiyiG\/TuSUVieRa2BdxGh8H8H8\/BM81qE4L6\/LCiE4OqiaSKi5BzTyVsQXN6HA8kIvHkJ8FtDjPMyj1pwCbXslapDylJ8i4Z4dVfpPVr1MIQst8qtJHfs2FG4b7MU3Cev+nzS3KlTp3ijJTJri4uLYu3atZHddfNuF7z52YN5ucy3bdt2mYH5+flyjZbQ+oqlpaWlUO3S8HFfX9+y5vImPOWJZ0dHR+M9rIngyrbku+S84wDRww3lcbN28O3fjFOoUuAdiqR5O2BuzipEyeh7uCEgqqJJP8uJVqZDylRHnXCVNsmKykBwQ3qruC0ko2JGIUuAd0iaZm2BuRmnUKUguKFIvj88vbCwsExw1R5vnimTSVoQ3IDOMmgKycgAUsAi4B0QpmFTYG4IKlCx6AWXepb9\/f3J0iD9Q8uDpqamRHt7uxFu12VBck1u0TIkCK6RG4IVQjIKhtKoIfA2whS0EJgHxVnYWNSCayp0hRSVAlkbX0hbvb29orOzU+gbX9hMmtp34SabS0JZRwJIRo7gHKuBtyM4j2pg7gHPoWrUgmu6HMeBaylV0MMtBWtmo0hG4M1LgN8aYpyXedSCq\/c6edHbW4Pg2jNs8DjoAAAbG0lEQVTzqYFk5EPPvi542zPzrQHmvgTt6kctuITKZLKSHdLySkvB\/dyhdeKm7mvKM4SWEwJIRryBAN68vBHj\/LyjE1w5jDw3N1dI23bSVGGDngUguJ4ALatDACyBeRYHb0+ADtXB3AGaR5XoBNeDVdOrQnB5XYBkBN68BPitIcZ5mUNweXl7WSPB\/ei6K5Lh5Fvvu86rLVQuJoBkVMwoZAnwDknTrC0wN+MUqlTUgps3S7mK73YhuKHC3qwdJCMzTqFKgXcokubtgLk5qxAlIbjvnxZEa2PVDwQ3RHjVuw0kI17\/gTcvb7IG5rzMoxRcdXOKPNyDg4NWB9CX7Tr0cMsmvLx9JCPw5iXAbw0xzss8SsGViOu48QXe4fI9IEhGfKzR2+JlLa0hxnm5Ry24vKj9raGH68\/QpgUkIxta\/mXB25+hbQtgbkvMrzwE148fa20ILituvN\/ixQ3ezLwxqsAPHILLz9zZIgTXGZ1TRXz7d8LmXAm8ndE5VwRzZ3ROFSG4TtiaU0kK7vquqwRt74hPuQSQjMrlq7cO3ry80cPl5w3B5WfubJEE9\/quq8TKdVdAcJ0pmleEAJizClESvENQtGsDzO14+ZaG4GYQpHW4fX19yV\/pYPjp6WmxceNGX95e9SG4XvisKyMZWSPzqgDeXvicKoO5EzbnShBcZ3T8FSG4vMyRjMCblwC\/NcQ4L3MILi9vL2sQXC981pWRjKyReVUAby98TpXB3Ambc6WoBbeOeynjHa5zrFtXRDKyRuZVAby98DlVBnMnbM6VILg120sZgusc69YVkYyskXlVAG8vfE6VwdwJm3OlKAW3znsp09F89MGyIOeYN66IZGSMKkhB8A6C0aoRMLfC5V04SsGV1Oq4lzIJ7juv\/0Lccezj3s5HA\/kEkIx4IwS8eXmTNTDnZR614JaBWu09j46Oiu7u7lQzUuzn5uaSv+\/YsUOMjY2Jtra2zMuiSVMQ3DK8lt4mkhEfayR\/XtbSGmKclzsEVwiRNsScJ5ZZLnr55ZfF0NCQGB8fT4rIn\/X1uxcvXhTDw8Oiq6srEWT5O6333bt3LwSX9xnItIZkxOsI8ObljS85\/LyjF1wS26NHj4qpqSnR3t6eeED2Pnt6ejJ7qGmuorZmZ2cbPdWHHnpIdHR0GLWh101rHz1c3gcEAgDevAT4rSHGeZlHLbihlwWRwNJH9lL13\/NcC8HlDXwTa0hGJpTClQHvcCxNWwJzU1JhykFwAy4L0nu0JKILCwu5w8Q2PWrq4Z798J+Lq5euFf\/pyq+LU6dOhYkCtJJKYHFxUaxduxZ0mAiANxNoxQyYl8t827ZtlxmYn58v12gJra9YWlpaCtFuyCFlF8GV72\/pXjBpKoRHw7WBb\/\/hWJq0BN4mlMKWAfOwPItai7qHK+GEmjRlO6RsI7Z0rdTDvXXPdeLHT78t7nnxxiLf4u+eBJCMPAFaVgdvS2ABioN5AIgWTUBwLWAVFdWHkPMmTZnOTFZtQnCLPBD270hGYXkWtQbeRYTC\/x3MwzPNaxGCG5C36bIgMklifP78+cJhZAhuQAdZNoVkZAnMszh4ewJ0qA7mDtA8qkQnuPqGE3nsNm\/evGy5kAnnrI0vZI+2t7dXbNq0SQwMDAi56YVst8geergmHghXBskoHEuTlsDbhFLYMmAelmdRa9EJrg4kb9IULe\/p7OwsYsj2dwguG+rEEJIRePMS4LeGGOdlHrXghl6HW7brpOC+MHFB7LtwU9nmom8fyYg3BMCblze+VPLzhuAGXIdbtvtIcElo5f\/Lthd7+xAA3ggAb17eEFx+3lELLuEOuQ63bPdBcMsmvLx9CAB48xLgt4YY52UeveAS7tOnT4u+vr5l5GdmZir1\/pYuDoLL+3AgGYE3LwF+a4hxXuZRCy7NHKb\/5KEFvOjtrUFw7Zn51EAy8qFnXxe87Zn51gBzX4J29aMW3KJTgQ4fPpycU1sVQYbg2gW3b2kkI1+CdvXB245XiNJgHoKieRtRC646nKweAE8bWPT394vVq1dbr8M1R29fEoJrz8ynBpKRDz37uuBtz8y3Bpj7ErSrH73gEi7Z033jjTcSkaUNKQYHBwtP+bFD7V8aguvP0KYFJCMbWv5lwdufoW0LYG5LzK88BPd9frJXS9stjo6OGh0a74fevjYE156ZTw0kIx969nXB256Zbw0w9yVoVx+C+\/6+xpOTk0mvduvWrcmMZXWI2Q5peaUhuOWxTWsZyQi8eQnwW0OM8zKPWnDVoeTp6WmxcePGhL78d\/p5amoKk6Z4Y7Iy1pCMeF0B3ry8yRqY8zKPXnBPnDgh7rzzzlTqmKXMG4xVs4ZkxOsR8OblDcHl5x214PLj9rOIIWU\/fra1IQC2xPzKg7cfP5faYO5Czb0OBNedHXtNCC4vciQj8OYlwG8NMc7LHILLy9vLGgTXC591ZSQja2ReFcDbC59TZTB3wuZcCYLrjI6\/IgSXlzmSEXjzEuC3hhjnZQ7B5eXtZQ2C64XPujKSkTUyrwrg7YXPqTKYO2FzrgTBdUbHXxGCy8scyQi8eQnwW0OM8zKH4GbwVo\/sW7NmjVDX6fK66ANrEFxe8khG4M1LgN8aYpyXOQSXl7eXNQiuFz7rykhG1si8KoC3Fz6nymDuhM25UtSCi\/NwneMmiopIRrxuBm9e3mQNzHmZRy24ZZyH+\/TTT4uRkZHEi6aHIDz00EOio6Oj8MAE2cN9\/Ja\/FPe8eCNvpERoDcmI1+ngzcsbgsvPO2rBJdzyXW2I83DpxKGhoSExPj6eeFL+LPdoTnMviS0dnGAizhBc3gcEAgDevAT4rSHGeZlHL7iEO9R5uNS7nZ2dFWNjY6KtrU3k9VylzZUrVyYev\/322616uDd1XyNuve863miJzBqSEa\/DwZuXN3q4\/LwhuO8zD3EeLgksffbu3Zv8X\/9ddS8J7ttvvy1oBvTw8LDo6uqC4PLHf65FCACvQ8CblzcEl583BDfgebh6j5Z6vAsLCw0BTnMvTdxyEdz1XVeJzx1axx8xEVmEAPA6G7x5eUNw+XlHLbihz8PlENwnr\/q8+IOff0u8u+KvxYEf\/SF\/xERkcXFxUaxduzaiO27urYI3P38wL5f5tm3bLjMwPz9frtESWl+xtLS05NsuCW7I83BthpTltbv0cJ\/6witi5bor0MP1DYCC+uhxlQxYax68eXmjh8vPO7oeruzV0nvWzs5OI+I0k5nE9NixY7nl9SFkk+U+EFwjFzSlEASAFzt48\/KG4PLzhuAaMDcVXJdlQS6C++xXXk+uGu9wDZznUQQC4AHPoSp4O0DzrALmngAtq0cruHNzc1aoNm\/eXNjDpQazNr6Qwtrb27usZ+0quO+8\/gtxx7GPW90DCtsRQDKy4+VbGrx9CdrXB3N7Zj41ohNcH1jNris3vqAeLgS3fG8gGZXPWLUA3ry8MaTMzxuCy8\/c2SIE1xmdU0UIgBM250rg7YzOuSKYO6NzqgjBdcLWnEpScF94+IL48dNvYz\/lkt2AZFQyYK158ObljR4uP28ILj9zZ4sQXGd0ThUhAE7YnCuBtzM654pg7ozOqSIE1wlbcyqpgvvCxAWx78JNzbmQSKwiGfE6Grx5eaOHy887esFV1+Vu2rRJDAwMCJrBTLOSp6amRHt7O79XMixKwaXhZJo4BcEt1zUQgHL56q2DNy9vCC4\/7+gFV90dipb0HD16NBHakydPFu6DzO0uCC4vcQgAePMS4LeGGOdlHrXgqr1b6tHSIQJ0eg\/tQiU3u6hSL1cKLoWI+jNvyMRjDcmI19fgzcsbPVx+3hDcgYFEYOVwck9PT3JMHgSXPxirZhECwOsR8OblDcHl5x214Kq7PK1fv17s2bNHTE9Pi40bN+aeZcvvpksW0cPlJQ8BAG9eAvzWEOO8zKMWXEKtHjw\/ODiY9Hbpve758+fF2NiYaGtr4\/VIjjUILq8rkIzAm5cAvzXEOC\/z6AWXF7efNQiuHz\/b2khGtsT8yoO3Hz+X2mDuQs29DgTXnR17TQguL3IkI\/DmJcBvDTHOyxyCq53wMzMzI86dOydmZ2cxpMwbi5WzhmTE6xLw5uVN1sCcl3n0givf1w4NDYldu3Yl73D1JUK8Lsm2hh4uryeQjMCblwC\/NcQ4L\/OoBTdtlykS3M7OzsovC3r8lr8UN3VfI2697zreiInIGpIRr7PBm5c3erj8vCG42jpcCC5\/EFbVIgSA1zPgzcsbgsvPO2rBJdy0nSO9r1WHlPVNMPjdkm5RHVJ+6guviJXrrhCfO7SuKpfXctcBAeB1KXjz8obg8vOOXnAJOe0q1dfXt4z+6OhosuNUlT4QXF5vQADAm5cAvzXEOC9zCC4vby9rquDSaUHvvP4Lccexj3u1icrZBJCMeKMDvHl5o4fLzxuCy8\/c2SIJ7pcfPSau+eMDyfF8EFxnlEYVIQBGmIIVAu9gKI0bAnNjVEEKRi246ixlmpmsfqp6eIEU3BceviDoXNx7XrwxSCCgkcsJIBnxRgV48\/JGD5efNwT3\/VnKzRBcmrA1MjKSeN3knTH1cO\/a\/12x+p7ppA6O6Cv3gYEAlMtXbx28eXlDcPl5Rym4qtDlIZeHGZThFjo0gWZGj4+PJ83Ln+mkoqwPBLcMT2S3CQEAb14C\/NYQ47zMoxRciThvSLlsN8jlSPJEItrxqqOjI3dmNAS3bK8sbx\/JCLx5CfBbQ4zzMo9acHlRL7dGAksf2miDPvrvaddGgvv3\/vl+8YXv\/Edx8Sd\/Jr77u58Rd\/+7e8RHVnWKq37n4aTK2\/\/+AXHxfz2X\/Nz2W7clE6zwcSOAZOTGzbUWeLuSc68H5u7sXGpGL7jqebg6QNpTeWpqSrS3t7uwza2j92ipx7uwsNAQ4CzB3dY9KH7\/wwfFc51D4tk\/nRLb9\/wP8YmuPxfvvXk6qfKrv3lP\/Opn74n3\/s\/PhfjlzeLnf3FJfOnz4dUd4sPXdoiP0P9Xd0CMC7yKZBQ87HMbBG9e3mQNzHmZRy24Fy9eFMPDw6Krq0vs3Lkz+bm3t1fInabkNo9luMRFcL\/3mavFx654T2y9+2Pixv+wV\/T\/qFP8cuOSeOOfvScO\/Nd\/IM587A\/EX2z8u+I3r\/ypWHvlT8VvXvlm8v+fzH5GnP\/Jb4t337xWvPvTa5P\/4wMCIAACIMBL4MmrPi\/m5+d5jQawtmJpaWnJtx39Ha4qgrQs6MiRI6Ud0ecypPzNT18j\/sXmleLXv7hWrPzHLyRrcX\/0yjviN1bfJ77waz8RO9fMiHXXXJlgub79SvGr0f8r1l34iPjouivER9d9BFtBWgYMvv1bAvMsDt6eAB2qg7kDNI8qUfdwdcFVh3XLXoerDyGbTJoiwb37T64UV\/7WbeLXP30kcTu91+353c3iN775Q9H2yduSf6M1ui9MXEiEFut03Z8OJCN3di41wduFml8dMPfjZ1s7asElWGpPUxXZkydPlnoIvcuyoH\/1968VB79xlWj7xFdF2yf+dUNw7\/infyiu\/zevJr9Tr\/fc7M8gtLZPQkp5JKMAEC2aAG8LWIGKgnkgkIbNRC+46ntcOqyABHhyclKsWbNGTE9Pi7x1sYaMM4vZbnxBznrxz5ZE+85L4ip7uLd\/\/Zdiy5\/+jqAzcumDXq2vZy7VRzIKw9G0FfA2JRWuHJiHY2nSUvSCawKpKmXSBHf+j1aI\/\/zq\/25cIsQ2nLeQjMKxNGkJvE0ohS0D5mF5FrUWteDWbS\/lsbs3iC\/9o\/fE3+l7PfErrbl997l\/K85tOC1uve+6Il\/j75YEkIwsgXkWB29PgA7VwdwBmkcVCG4T91K29dvhr28Qv7\/+5+Lau\/4qqXr+m7+HzS1sIVqURzKygBWgKHgHgGjZBJhbAvMsHqXgVmEvZRe\/nfzOBrH56g8El4aTb3jGe3WUy6VEUQfJiNfN4M3Lm6yBOS\/zKAVXIm7mXsoubibB\/eTP3hFX3\/pdcfVtfyJeu3tDY3ayS3uok08AyYg3QsCblzcEl5931ILLj9vP4lvHN4i\/+S8iWW\/7\/95YSBpb880f+jWK2pkEIAC8wQHevLwhuPy8oxVcGlZ+9NFHG0t\/ZG93bm4u8YLJ+bTc7iLB\/eXrlza3oMlSGE4u1wMQgHL56q2DNy9vCC4\/7ygFlza42LNnT0Ns5VpcWntL+ydL8e3p6ck9Lo\/bXdJZNJRMH7nZBfd1xGIPAsDrafDm5Q3B5ecdneDqG10Qcnli0MTEhOjs7Ey8oJ9Xy++ayy3W1VlVYOdyDRAAF2rudcDbnZ1rTTB3JedWr6453PnwgrSJUnqPl1CWvZeyi7vq6iyXe61CHSQjXi+ANy9v9HD5edc1hwcV3LTeLASXPxirZhECwOsR8OblDcHl5x2d4MohZTr3loaP9fe30gW0p\/L58+dLO57PxdV1dZbLvVahDgSA1wvgzcsbgsvPu6453LmHS4jVHu3i4qLo7+8X6vtb6t329fWJmZmZxjtdftdcbrGuzqoCO5drgAC4UHOvA97u7FxrgrkrObd6dc3hXoJLqOSpQPSzXAIke7tnzpwp\/aQgF3fV1Vku91qFOkhGvF4Ab17e6OHy865rDvcWXH7U\/hbr6iz\/O29OCxAAXu7gzcsbgsvPu645HILLHyvRWYQA8LocvHl5Q3D5eUNw+Zk7W6yrs5xvuMkVIQC8DgBvXt4QXH7edc3h6OHyx0p0FiEAvC4Hb17eEFx+3hBcfubOFuvqLOcbbnJFCACvA8CblzcEl593XXM4erj8sRKdRQgAr8vBm5c3BJefNwSXn7mzxbo6y\/mGm1wRAsDrAPDm5Q3B5edd1xyOHi5\/rERnEQLA63Lw5uUNweXnDcHlZ+5ssa7Ocr7hJleEAPA6ALx5eUNw+XnXNYdXtodL20aOjIwknjQ9xJ52vero6Cg8e7euzuIP6zAWIQBhOJq2At6mpMKVA\/NwLE1aqmsOr6Tg0rm6Q0NDYnx8PGEvf964cWOmL+QWkybiXFdnmQRiFcsgGfF6Bbx5eaOHy8+7rjm8koKrH\/OX13OV5\/KuXLky8frtt9+OHi5\/\/OdarOvDUTGMxpcD3saoghUE82AojRqqK+9KCi4JLH327t2b\/F\/\/XfUICe7bb78t1qxZI4aHh0VXVxcE1yhk+QrV9eHgIxTWEniH5WnSGpibUApXpq68Kyu46rtY6vEuLCw0BDjNbfKEIlPBDed6tAQCIAACIMBNYH5+ntukt70oBdebGhoAARAAARAAAUsCTRdcdTby5s2bxdTUlHjiiSeMh5Tl\/dr0cC0ZoTgIgAAIgAAIeBNouuCm3YE+hGyy3AeC6x0LaAAEQAAEQKBEApUUXJdlQRDcEqMETYMACIAACHgTqKTg0l1lbXwhhbW3t1d0dnY2AEBwvWMBDYAACIAACJRIoLKCW+I9o2kQAAEQAAEQYCcAwWVHDoMgAAIgAAIxEoDgxuh13DMIgAAIgAA7gegEV+65TKRnZmaWvQdmp9\/iBuW2m3Nzc407NdnrusWxlHJ7aTP5Vf5yyV17e3sp9mNrNI336dOnRV9fXwMF7X43PT0t8vaAj42bzf3q+WPHjh1ibGxMtLW1Jc3UMb6jElx6IOhBobW+Z8+ebfyMJGTzGJiXpdnmBw8eFI888ogAY3NutiWzDu5Qt0TN2x7V1l7s5bN4m+yIFzs70\/vXJ8HK3+lLTNqWv3WJ76gEV3VK1mxn04BAuWIC9AXnyJEjy76VFtdCCVMCeQd3yL9RcqLZ\/PjyY0o1u1zRQSkm+wX4X0W8LaiH2lD+HhgYSMS3TvEdjeBmfWMy2Xs53hD3u3P91Ce\/1lBbJ5B3cIe6lp2GNPXfQdOeQB5vLEu052lbQ80ni4uLy45trUt8Rye46vpdfCO1DXm78ur7cqqpv4Oxaw2lswikJXu9R0tisXv3brFv3z68U\/QMpTTemK\/gCbWguuTb09OTnAZX1\/iG4HZ0FB7nV24otWbr+juXtHcwrXnn\/HcFweVlnsW7v79fTExMNIY41d95r7C1rEnedFdy0hQEt+I+xpBy8x2kTlrDJKpw\/sgSgKGhITE+Pp70aOsy5BaOSnktmQ4f12UiT3mk\/FtOE1tqta6vTKLp4ZKT1CFkTJryfxhsW8AkKltiZuWzhjjVIWRMmjJjaVLKRnDVc71N2kaZDwjkjYrpr0jqEt9RCS6WBfE9zvosWf0dDN+VtL6lLAHAsqByfJ\/GWx+9od\/37NmDdbgeLqD4PX\/+fOYqhzrGd1SCK3u5k5OTSRhg4wuPp8Ggqj6RZHBwsLGGzqA6ihgSyBLcOm4MYHjLTS2WxVvf+AL5xd1NaZPQqDV1A5c6xnd0guseAqgJAiAAAiAAAu4EILju7FATBEAABEAABIwJQHCNUaEgCIAACIAACLgTgOC6s0NNEAABEAABEDAmAME1RoWCIAACIAACIOBOAILrzg41QQAEQAAEQMCYAATXGBUKgoAfAVqcT9v90drCrM\/DDz8snnrqqcYOUX4WzWubbkpSlw0GzO8cJUGAjwAEl481LIHAMgJV2frP9mADnPuKQAYBNwIQXDduqAUC3gSqIri2Amor0N6g0AAItAgBCG6LOBK3UT8CaYKrbsq+atWq5JDtnTt3iqmpqWQoes2aNcl2gS+99JIYGRlJblo\/9lDf8Wh0dDTzRCx9C05qT+6kdOLEiaR9aZMOQZAfW5Gun3dwxSAQngAENzxTtAgCRgRMBZcaI8GlE5bkGcNSRPU9qkkIjx492ihftId12glO+nWpB3+3tbUl94aTn4xcjEIgsIwABBcBAQJNImAquPLQ7TShU\/f1pZ7w8PCw6O3tTc5klZ88cdTF1PTcYhz316SggdlaE4Dg1tp9uPg6EzAV3L179zYEVBdPVXC3b9+eDEHPzc1dhkXd9F39Y1bvta+vLymWNpxM\/542FF1nX+DaQYCDAASXgzJsgEAKgdCCu2XLlmTZ0cTExLIebh78NMGV5dV3ubrwQnAR0iBgTwCCa88MNUAgCIHQgit7uOoQdNGF5gmuLrxdXV2NyVcYUi4ii7+DwOUEILiIChBoEoHQgtvd3S1IQB999NFlB5\/nHeStb2SR9g43TVwxaapJQQOztSYAwa21+3DxdSZQhuASDxJduWSIfteXDanMpMCqE63SDv\/WD1OvyhriOvsf1x4fAQhufD7HHYPAMgK2a2qx8QUCCATcCEBw3bihFgi0DAFbAbUV6JYBhRsBAU8CEFxPgKgOAq1AAIcXtIIXcQ9VJwDBrbqHcH0gAAIgAAItQeD\/A6le19FlVbIyAAAAAElFTkSuQmCC","height":286,"width":476}}
%---
%[output:9f4e2fbc]
%   data: {"dataType":"text","outputData":{"text":"Closed-loop eigenvalues:\n","truncated":false}}
%---
%[output:2cb5c6d2]
%   data: {"dataType":"text","outputData":{"text":"   1.0e+03 *\n\n  -1.0655 + 0.0000i\n  -0.0091 + 0.0027i\n  -0.0091 - 0.0027i\n  -0.0046 + 0.0000i\n\n","truncated":false}}
%---
