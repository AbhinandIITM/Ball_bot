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

% --- 3. Derived Parameters ---
l = h + r; % Height of ball CoM from wheel center


M1 = [(I_b/(r^2)) + m_b, m_b*l;
      m_b*l,            I_r + m_b*(l^2)];

% M2 = [[1,  1+(g*m_b)],
%       [1-(m_b*g), 1 - (M_r*d*g)-(m_b*g*l)]]
M2 = [1,             (g*m_b);
      (m_b*g), 1 + (M_r*d*g) + (m_b*g*l)];

% M3 = [[m_b*l],
%       [-m_b]]
M3 = [-m_b;
      -m_b*l];

% --- 5. Calculate A and B ---
% This is the numerically preferred method for A = inv(M1) * M2
A1 = M1 \ M2;

% This is the numerically preferred method for B = inv(M1) * M3
B1 = M1 \ M3;
%%
A = [0 1 0 0;... %[output:group:1633ba1f] %[output:60988226]
    A1(1,1) 0 0 0;... %[output:60988226]
    0 0 0 1;... %[output:60988226]
   A1(2,1) 0 A1(2,2) 0] %[output:group:1633ba1f] %[output:60988226]
B = -[0 ; B1(1); 0 ; B1(2)] %[output:06e878d2]
C = [0 0 1 0] %[output:717ff987]
D = [0] %[output:18f3fd5a]
%%
[num,den] = ss2tf(A,B,C,D) %[output:052d9a1f] %[output:17d34675]
%%
G = tf(num(1,:), den) %[output:7d01220c]
%%
sys = ss(A, B, C, D);
%%
rank(ctrb(A,B)) %[output:4777e82b]
rank(obsv(A,C)) %[output:4dc76433]
%%
x0 = [0.005, 0, deg2rad(0.2),0]' %[output:6e1401a3]
%%
OS = 10;      % percent overshoot
ts = 2;       % desired 2% settling time (seconds)

% correct 'a' and zeta
a = -log(OS/100);            % a = -ln(OS/100)  (natural log)
zeta = a / sqrt(a^2 + pi^2); % damping ratio

% natural frequency for 2% settling time approximation: ts ≈ 4/(zeta*wn)
wn = 4 / (zeta * ts);

roots([1 ,2*zeta*wn , wn^2]) %[output:50e88050]


poles = [ -2 + 2.73i, -2 - 2.73i,  -60, -40 ];

K = place(A,B,poles) %[output:02d03c09]
% K = a 1x4 row vector

% Closed-loop A matrix
Acl = A - B*K;

% Check poles
disp('Closed-loop poles:') %[output:1f9316b1]

%%
sys_cl = ss(A - B*K, B, C-D, D);

figure; %[output:4c5896db]
step(sys_cl); %[output:4c5896db]
title('Closed-loop Step Response (Pole Placement)'); %[output:4c5896db]
%%
x_inf = -(inv(A-B*K))*B %[output:47d5323b]
u_inf = -K*x_inf + 1 %[output:8ef96ea5]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:60988226]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["263.6965","0","0","0"],["0","0","0","1.0000"],["-0.4806","0","41.5546","0"]]}}
%---
%[output:06e878d2]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["0.7141"],["0"],["0.0006"]]}}
%---
%[output:717ff987]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":1,"type":"double","value":[["0","0","1","0"]]}}
%---
%[output:18f3fd5a]
%   data: {"dataType":"textualVariable","outputData":{"name":"D","value":"0"}}
%---
%[output:052d9a1f]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":1,"type":"double","value":[["0","0","0.0006","-0.0000","-0.4966"]]}}
%---
%[output:17d34675]
%   data: {"dataType":"matrix","outputData":{"columns":5,"exponent":"4","name":"den","rows":1,"type":"double","value":[["0.0001","-0.0000","-0.0305","0.0000","1.0958"]]}}
%---
%[output:7d01220c]
%   data: {"dataType":"text","outputData":{"text":"\nG =\n \n            0.0005817 s^2 - 2.067e-18 s - 0.4966\n  --------------------------------------------------------\n  s^4 - 4.441e-15 s^3 - 305.3 s^2 + 2.274e-13 s + 1.096e04\n \nContinuous-time transfer function.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 78 117 109 101 114 97 116 111 114 58 32 123 91 48 32 48 32 53 46 56 49 54 55 101 45 48 52 32 45 50 46 48 54 54 53 101 45 49 56 32 45 48 46 52 57 54 54 93 125 10 32 32 32 32 32 68 101 110 111 109 105 110 97 116 111 114 58 32 123 91 49 32 45 52 46 52 52 48 57 101 45 49 53 32 45 51 48 53 46 50 53 49 49 32 50 46 50 55 51 55 101 45 49 51 32 49 46 48 57 53 56 101 43 48 52 93 125 10 32 32 32 32 32 32 32 32 86 97 114 105 97 98 108 101 58 32 39 115 39 10 32 32 32 32 32 32 32 32 32 73 79 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 48 10 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 39 39 125 10 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 39 39 125 10 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:4777e82b]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:4dc76433]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"4"}}
%---
%[output:6e1401a3]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x0","rows":4,"type":"double","value":[["0.0050"],["0"],["0.0035"],["0"]]}}
%---
%[output:50e88050]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":2,"type":"complex","value":[["-2.0000 + 2.7288i"],["-2.0000 - 2.7288i"]]}}
%---
%[output:02d03c09]
%   data: {"dataType":"matrix","outputData":{"columns":4,"exponent":"5","name":"K","rows":1,"type":"double","value":[["0.0462","0.0017","-3.0913","-0.3189"]]}}
%---
%[output:1f9316b1]
%   data: {"dataType":"text","outputData":{"text":"Closed-loop poles:\n","truncated":false}}
%---
%[output:4c5896db]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAcoAAAETCAYAAAC2rUY7AAAAAXNSR0IArs4c6QAAIABJREFUeF7tnX1sXled58\/s7Cz1FGYTZwNMYiBV66CuVsqUbFU3rRYQgyIxcroqq3Fi\/kBeE1miJUFN4peEF1FB4jjNSGnLoCgYC2mxk0UN21paBnVYRJaGRhUvYf+oNqFMahxPUaduV8Ak7O4sq99Nz5Pj4\/t27svz3HPv50qW\/Tw+99xzPr\/fPd\/7O2\/3D37\/+9\/\/XnFAAAIQgAAEIBBK4A8QSjwDAhCAAAQgEE0AocQ7IAABCEAAAjEEEErcAwIQgAAEIIBQ4gMQgAAEIACBbASIKLNx4ywIQAACEGgIAYSyIYammhCAAAQgkI0AQpmNG2dBAAIQgEBDCCCUDTE01YQABCAAgWwEEMps3DgLAhCAAAQaQgChbIihqSYEIAABCGQjgFBm48ZZEIAABCDQEAIIZUMMTTUhAAEIQCAbAYQyGzfnsy5fvqyGhobU0tJS69wjR46ogYGB1uejR4+qkydPqtnZWdXX1+d8jaJOeP7559Xg4KAaGRlRY2NjodlWpaxm4XS59Xf9\/f1qcnJSdXV1FYUmVT6aTVjiTpUpVcE7kEjfF5\/61KeCeyGK3YYNG9TMzIzq7e1NVcrl5WU1PDwcpJ2enlbd3d2pzjMTpSlLEddxLliFThBGct9lZVyhqsQWBaFsg6XOnDmjJiYmQq9kNpxVER8fhdIWySixFMbvf\/\/7S30QiRNKKVfcA0gb3LFSl7Ab2jh2LmJZhIClKcu6detyC3KlDJJQmGvXrqkvfelL6uMf\/3jw0GI\/6PhUF5eyIpQutDKkNSNJM1I0xVN\/j1BmAPzmKZqdjtI1d\/m3RCI9PT1qfHxczc\/Plx6xR9lRl2n9+vW1fwJPY0nNQx4Wdc9FGDtpnLXt7F6YqOsUKZTmfWuX5cMf\/nBjhFIzffXVV1vRveaxsLBQa59GKNPc0TnSaEEMu8GlUdi0aVOr+zWskbC7bLds2bLCIc0bVxfTjljsiNbu2jWjMcl\/9+7d6uGHH3buetU30sWLF4OihEUASfXRDJ588kl16tQppfNK6o7W54V1bYYx0umknLoRlr9tvlnKEyWUYQ2Nzcwuf1iXvWZhioG2WVgd5DvbB0wfMfMRwdq\/f39riMDkHlcWuYbN2WZp30Zh90YUO53WLHfaOuluQdfypSmLcLe7eNMMs9g9IPY9a\/\/fbD90uR555BH13e9+N7hH9L0mjPUQj+1LcXnqMm\/dulU98MAD6hOf+MSKe1hHzvp+NHtGqvKAn6OZTjwVoUxElC+BixPZaaO6E\/VNYUZJdin1jRXV7asbwKhrmDdCGAG7rGGNgz4v6VqmoMZ1d8WJZVzDEyWUn\/vc59Sjjz4aRJnmYTbwWcoTZXNdRp2\/XFMaWbPxke90A7e4uLhqXNt8AAlrvHQ90tRBX0f4hJXDFN3XXnsttixRvhgnlmGc0kaUUXax6yR1EKGUcWrzgSiMk+3nacpiR5RRNjV9POqe02IZ9X99Tyd17Zv10Ock5Rl3\/wrT0dHR4OE5TCjjgoF8rWd1zkYoS7SF2UAnRURSDPPGlAZG39jm06ZOI9\/t2bMnSPOjH\/0odKJDWARjPjnKRJfHH388mECkr2GW2WUyT1hEp28gudFMUYqqj0QzYfnERYum+cIeCnRDEWYL3XiYT972TZ+lPEkNmf0QE8ZD\/EWEUKIDecoPm5RkRqM6T\/u7973vfS2BsyNRafTku82bN7eEMqrrWkcqUWVJw9K0lbaH7btx7MKipqQ6aaG8dOlSMEEtztZRQhnWROiyJI1R2n5n3te2zeQ60pMyNTW1YojA7kaWnha5Z8MeqMJ8S7cT5rCDnaf5IGQ\/2NoPdmbXq5Q5zZyGEpvZtmSNUJaIOY9Q6sZL3+h61p49zvXss8+umChkPrnGPSWK8+ub0m6sTMc3bzIzQvz+97\/fmqGrb347H1OoT5w4ob74xS8GWZgz5Oz66EYgrMsvzdheWPQoeZkNlM47bpKVFq+47vCo8kQ19nZXdJwoSIO3Y8eOVVGQ2dBHjcOZ3ZQycSlsBrOZRncfho09aZuGRYxhohN2O4U9cIU9xJkPi3Y+ZmQa1TCH1Un7m32fmPlHPRBG2ccsS5QNovww6r7W5bG74s1yav85e\/ZscO\/ZQmvaz2Sk7Wv3XEje9sOH6dN23XS0HCWUdZ7RjVCWKJSSddjYir6k3IhRY5RphVILaNh4zYMPPhjaXSbXr6tQarb2pAtTdNollGFdzmFjTWEuaKZzbXSzCqX5EBMV8UWV5eWXX46c2Z1FKNN0s0eNxcv39thhHqGMK0uUmGiWZpevGb3bD4xFCKWZZ1ahNHsMEMqbdyZCWbJQRs16NccMwma9pul6DVvjaEZox44dUwcOHFD2E6BZZbMrV\/KrYter2YUb1wVp11PXzR7bsQUs7kk4qUs5rDxhUWiYveMeoqLcMixqkkghrItZvnPpek0jlGa5wgQ5bVSR1PUaJ05h95QZiYUJUljXa9Ktn2Z+gS0mWpC1iNvliut6Ff\/9yle+or761a\/Gzs62\/Tosqo3qFYriag\/JiMDb+WrRj+uBilp3ncS66v9HKNtgoazrKJMm88RN5tA3atS1kyYOCBaXMcpOT+ZJ091pprHHTW03sB9ewtwkqtGJamD190mTeXR3mIwbSbdp1LXNsUU7TZbJPHFCmVQWUwCiWNrfp53ME1b\/oibzSN6udjTLY4uJFuQ4f0la8ysPPmF21w8hel5B2BinHtawu6eTrqknjsVFlPaEKLuNSbt0pw1NbuGXQCgLRxqeYZop43HjYXpHnzTLQ+yn+rhp9FLaKi4Pkafqp59+ujUjNc1kKLue9pigXc+o2ZDmtbRNXMoTJZRmdBEWcYgt7DLH+U3Uso6wmaZZl1KY0UOSD0d1y0btMuWyPCTsrspaJ3OWc5xfZYkoZSjEFHGxc9g4cdIuUnH\/d40odZQXl2eaiFLqZjK3hTvNPdqm5rbwyyCUhSMlwzwE0jROefJ3Pbdq5YmLZrJs0+bKo8j0YRsOFJk\/eZVPgA0HymfMFSCwikDVhKlq5amTUEpdhK9EOnXfK7Sutzpb2NXVstSr0gSqJkxVK0\/dhLIpDW2lb7ochWvKgw5drzmchFMhAAEIQKD+BBDK+tuYGkIAAhCAQA4CCGUOeJwKAQhAAAL1J4BQ1t\/G1BACEIAABHIQQChzwONUCEAAAhCoPwGEsv42poYQgAAEIJCDAEKZAx6nQgACEIBA\/QkglPW3MTWEAAQgAIEcBBDKHPA4FQIQgAAE6k8Aoay\/jakhBCAAAQjkIIBQ5oDHqRCAAAQgUH8CCGX9bUwNIQABCEAgBwGEMgc8ToUABCAAgfoTQCjrb2NqCAEIQAACOQgglDngcSoEIAABCNSfAEJZfxtTQwhAAAIQyEEAocwBj1MhAAEIQKD+BBDK+tuYGkIAAhCAQA4CCGUOeJwKAQhAAAL1J4BQ1t\/G1BACEIAABHIQQChzwONUCEAAAhCoPwGEsv42poYQgAAEIJCDAEKZAx6nQgACEIBA\/QkglPW3MTWEAAQgAIEcBBDKHPA4FQIQgAAE6k8Aoay\/jakhBCAAAQjkIIBQ5oDXjlOff\/55NTg4GFxqdnZW9fX1teOyXAMCEIAABN4kgFBW2BWWl5fVqVOn1J49e9S1a9daf3d1dVW41BQNAhCAQL0IIJRtsOfly5fV4cOH1fHjx1V3d3dwRRHB4eFhdfHiRbVlyxY1PT3d+p8ukpx39uxZNT8\/r5aWlogo22ArLgEBCEDAJoBQluwTInZDQ0Nq\/fr1K8Tw6NGjwZXHxsaU+bdZHDn361\/\/ujp06FAQUX7hC19Qn\/\/851cJaslVIHsIQAACjSaAUJZo\/jNnzqiJiQk1MjKiZKxRR406mhSRlDFHM+K8dOlSMCa5YcMG9eijj6rr16+rj3zkI4FQPv7442r37t0IZYk2I2sIQAACRJQF+oCI38svv6wGBgaCXOXz3NycmpycVDKOKN2qmzdvDn5L1KiFUoRxdHRUTU1Nqd7e3kAozc+6iCKoOoqU7\/R4JWOUBRqRrCAAAQgkECCizOEiEuWNj4+rhYUF9aEPfSgQPC2SZrYioLZQmmOWIoj79u1TBw8eDITTPpdZrzmMxKkQgAAEchJAKHMClNOjxhh11nmEsoDikQUEIAABCOQggFDmgJcnokzT9ZqjaJwKAQhAAAIFEUAoc4BMGqOMiijtrtaw5SM5isWpEIAABCBQIAGEskCYUVnZXa92d21S120bisglIAABCEAgggBC2QbXCBPKNBsOtKFoXAICEIAABBIIIJS4CAQgAAEIQCCGAEKJe0AAAhCAAAQQSnwAAhCAAAQgkI0AEWU2bpwFAQhAAAINIYBQOhp6cXHR8QySQwACEICASaCnp8crIAilg7lEJA8cOKAuXLjgcBZJIQABCEDAJHDPPfeoY8eOKV8EE6F08F9Z5iH7roqBN27c6HBmPZPKA8OJEyfgoVTw8ASLm34Oj5X3PDxW+8a5c+cQyjpKgxZKnwxcph3gcZMuLFZ6GjzgEdX2+OgbRJQOSuKjgR2q55xUuqKfeuoptXfvXudz63YCLFZaFB7wQCjr1sqlrA9CmRIUySAAAQhEEPCxHSWidHBnHw3sUD2SQgACECidgI\/tKELp4BY+GtiheiSFAAQgUDoBH9tRhNLBLXw0sEP1SAoBCECgdAI+tqMIpYNb+Ghgh+qRFAIQgEDpBHxsRxFKpdSZM2fUxMRE4CBHjhxRAwMDoc7io4FL93ouAAEIQMCBgI\/taOOF8vLly2p0dFRNTU0FptZ\/9\/b2rjK9jwZ28F+SQgACECidgI\/taOOFUqLJ8+fPq8nJSdXV1aWOHj2qNm3aFBpV+mjg0r2eC0AAAhBwIOBjO9p4oRRhlGNsbCz4bX827R9l4IXl6+rd3bc4uApJIQABCDSTAELpod3tCFIizCtXrrSE0xbKnZ\/Yo\/7i0NeViKMcz730xopaa8Hcdfc71a67\/xQB9dAnKDIEIFA8Af3mJfkte2b7tBUoEaXV1ZpGKH9z\/wH1z\/\/hfwaeNPD+f6P6+vrUu9beiCife+n1G5Hpd64Ev0U4RTTHtt9WvOeRIwQgAAFPCMhLA+RHHwilJ4YLBC1D16v59hB5TUzUq2Ik6px74e8D0RTBfHLnner+O9Z4RIeiQgACECiGgESS8qPfpIJQFsO1LbnYEWQZk3lEMB+aezHoph3bvonosi2W5SIQgEAVCTBGWUWrJJSpnctD5l54JRBMoksPHYUiQwAChRBAKAvB2P5M2rnhgBldPvPJu+iKbb+5uSIEINBBAghlB+G349JFGlgiS4kwv7zrzmCyDwcEIACBJhAosh1tF6\/Gz3p1AV20gXVXrAilCCYHBCAAgboTKLodbQcvhNKBchkG1mLJJB8HQ5AUAhDwlkAZ7WjZMBBKB8JlGfgHP39D7fjrnwRdsESWDgYhKQQg4B2BstrRMkEglA50yzQwYulgCJJCAALeEiizHS0LCkLpQLZsA+tuWCb4OBiFpBCAgFcEym5Hy4CBUDpQbYeBj37n74KdfFg64mAYkkIAAt4QaEc7WjQMhNKBaLsMrJeOIJYOxiEpBCDgBYF2taNFwkAoHWi208D9X\/5JsOUdYulgIJJCAAKVJ9DOdrQoGAilA8l2G1jE8pevX1c\/\/cy9DqUkKQQgAIHqEmh3O1oECYTSgWK7DSzb3f3ZF3+o7rt9jZp\/6C6HkpIUAhCAQDUJtLsdLYICQulAsRMG1stG2JDAwVAkhQAEKkugE+1oXhgIpQPBThmYmbAORiIpBCBQaQKdakfzQEEoHeh10sB6co+MV8prujggAAEI+Eigk+1oVl4IpQO5ThtYj1eyzZ2D0UgKAQhUikCn29EsMBBKB2qdNjDjlQ7GIikEIFBJAp1uR7NAQSgdqFXBwIxXOhiMpBCAQOUIVKEddYWCUDoQq4qBWV\/pYDSSQgAClSJQlXbUBQpC6UCrKgZmfaWD0UgKAQhUikBV2lEXKAilA60qGViPV7LFnYMBSQoBCHScQJXa0bQwEMq0pJRSVTMwS0YcjEdSCECgEgSq1o6mgYJQpqH0ZpoqGliWjLxr7S1scedgR5JCAAKdI1DFdjSJBkKZRMj4fxUNzJIRBwOSFAIQ6DiBKrajSVAQyiRCFRdKKR5LRhyMSFIIQKCjBBDKjuIv\/+JVNrCMV8rBW0bK9wOuAAEIZCdQ5XY0qlaNjyivXbumxsfH1fz8fIvRyMiIGhsbW8WsygamCzb7jcuZEIBA+whUuR1FKCMILC8vq3379qmDBw+q3t7eWG+puoHpgm3fzc6VIACBbASq3o6G1Sp3RHnmzBk1MTER5D07O6tefvlldf78eTU5Oam6urqykWzjWZcvX1aHDx9Wx48fV93d3V4LpRSeLtg2Og+XggAEnAk0TiiPHj2qlpaW1OjoqHr44YeD7sotW7YEXZkbNmwI7b50plryCWI0qcf09HRqoZQHgp6enqBk+nfJxUydPV2wqVGREAIQaCOBxcXF4Grye3BwUJ07d65y7WcUjswRpXRZDg8PB2K4efPm1t99fX3Bwvy04tNGO4VeyoyIJYEIfZRo6ichM6O9e\/cq+anSobtgeXdllaxCWSDQbAInTpxQ8qMPhNIjodRRse4qtj+brq2F8tixY2rjxo2tiLJqUaUUjC7YZjdK1B4CVSMgkaT8XLhwIRDMRgilGEGiMRmPNLtedXS5c+dONTAwUClbyXjk0NBQ0F0sh3ShSgRsHpJG6jM1NbVqco9Pfeu6C1Ze8rzr7ndWyg4UBgIQaC4Bn9pRbaXMXa86g7DuyCNHjlROJNO6ZdzkHt8M\/NDci2ruhVcUXbBprU86CECgbAK+taPCI7dQlg21zPz1Gspt27YFwq4\/R01E8tHA7AVbpgeRNwQg4ErAx3a00UIpBrY3HOjv749c2uKjgXkdl+ttTHoIQKBMAj62o05CqWe6Xrx4MZFj3OzRxJMrmsBHAwtK\/Tqu5b\/6YEXJUiwIQKApBHxsR52E0jakTOY5ffr0iuUU5rIRe6KM747go4E18+5HvhdM6pHJPRwQgAAEOkXAx3Y0s1DGCaJP6yhdnMVHA+v6yaQemdzzzCfvUvffscal2qSFAAQgUBgBH9tRhNLB\/D4a2KyedMH+8vXrwSxYDghAAAKdIOBjO5pZKAVwXNdrFddR5nUKHw1s1nlh+bqSWbBj2zepse235cXB+RCAAAScCfjYjuYSSiEUto4ybCG\/M80KnuCjgW2MvGGkgo5FkSDQIAI+tqO5hbJB9m09FPi09VKYfdjerkleS10hUC0CCGW17FF4aXw0cBgEtrcr3DXIEAIQSEnAx3Y0c0SZtKaSdZQpvaZDydjerkPguSwEGk6gUUIZZWsR0H379qmDBw+u2lTcd\/\/w0cBxzGViz323r2Ftpe+OSfkh4BEBH9vRzBFlnF0ExNzcXORWcB7ZdEVRfTRwHGvWVvrqiZQbAv4S8LEdLU0ofXlxs4u7+WjgpPqxtjKJEP+HAASKJOBjO1qKUMa9\/LhI4O3Oy0cDJzFibWUSIf4PAQgUScDHdjSzUMZN5pHXVM3MzDBGWaR3lZgXaytLhEvWEICA90NYmYWyibb38UkorZ1YW5mWFOkgAIE8BHxsRzMLJZui9+Txlcqdy9rKypmEAkGglgQQyjfNyttD\/PRv1lb6aTdKDQGfCDRCKGUj9ImJiUS7jIyMqLGxscR0PiXw0cCufGVt5bvW3qLmH7rL9VTSQwACEEgk4GM7WkrXayIpTxP4aGBX1LoLlvdWupIjPQQgkIaAj+1oZqFMA6RuaXw0cBYbyMSe5156Qy3\/1QeznM45EIAABCIJ+NiOOgmlOYFn8+bNanh4WF28eDEUCHu9+n2ndD\/yPbXr7neyvZ3fZqT0EKgcgdoLZeWIt7lAPho4KyK2t8tKjvMgAIE4Aj62o04RZdPN76OB89iM7e3y0ONcCEAgjICP7ShC6eDLPhrYoXqrkrK9XR56nAsBCDRSKJPeQWlCYYyyHjcJ29vVw47UAgJVIeBjwNGoiPLatWtqfHxc7dq1S\/X19bX8xlwbeuTIETUwMBDqUz4auIibg+3tiqBIHhCAgBDwsR1tjFBqkZyfn1ezs7Mtobx8+bIaHR1VU1NTgRfrv3t7e1d5tY8GLuLW1Gsrx7ZvUmPbbysiS\/KAAAQaSsDHdjS3UIbt1GMKURV8QcRwaGhIbd26VS0sLAQ7BumIUsp\/\/vz51kum5RVhmzZtCo0qfTRwUfzpgi2KJPlAoNkEfGxHcwmliMzp06fV9PS06u7uDqyvxzF37twZ2YXZbje5evVqcMmurq5g7acplCKMcujt9uzPZlm1geVBoKfnxqbo+ne769SJ69EF2wnqXBMC9SCwuLgYVER+Dw4OqnPnznnTfmYWSh\/fHhJWZjuCFPG\/cuVK6D61WihNt927d6+SnyYcvGGkCVamjhAoh8CJEyeU\/OgDoXz+eSUCZEaa5aB3y7UooTx27JjauHFjK6JsUlTJG0bcfI7UEIDADQISScrPhQsXAsFshFBKxSXC2r9\/v5qZmVF68kunu15FoE+ePBkYpr+\/vzX2KJ+jhFL+59L16pOBy7hJecNIGVTJEwLNINCoMUrXNZXf+ta3Ou4FYUJpd7UymSfZTHTBJjMiBQQgEE6gUULpoxOECSXLQ7JZki7YbNw4CwJNJ9A4odTLLpaWllbZvoo780RNQGLDgWy3Ll2w2bhxFgSaTKBRQqkX8G\/YsCF0hmgdHcFHA5dpB7pgy6RL3hCoJwEf29FSlofU07x+br1Uti3ogi2bMPlDoF4EGiWUUfum1sukK2vjo4HbYQ82ImgHZa4BgXoQ8LEdzRxRismkwlVcL1mWO\/lo4LJYmPmyF2w7KHMNCNSDgI\/taC6h9G0yT14389HAeeuc9nz2gk1LinQQaDYBH9vRzEKpu163bdtWmT1dy3Y\/Hw1cNhMzf7pg20mba0HATwI+tqOZhTJur1c\/zZdcah8NnFyr4lIsLF9XsmSE13EVx5ScIFA3Aj62o5mFksk8N94ewrGSwNwLryiZCfvMJ+9S99+xBjwQgAAEVhBolFBKzWWM8vDhw+r48eOt12zV2Sd8NHAn7CFdsL98\/br66Wfu7cTluSYEIFBhAj62o5kjyqS9Xqu4M09e3\/HRwHnrnPX87ke+p+67fY2af+iurFlwHgQgUEMCPrajmYWyhvZLrJKPBk6sVEkJ9JIRumBLAky2EPCUgI\/taGlC+dJLL6m1a9fWqkvWRwN38l5iyUgn6XNtCFSTgI\/taKFCqSf4zM\/PK7peq+mk7S4V45XtJs71IFBtAo0VSvPtG2KikZGRWm6U7qOBO33L6CUjjFd22hJcHwLVIOBjO5o5ogybzFNXgdTu5aOBq3BryHjlw6dfVE\/uvJMlI1UwCGWAQAcJ+NiOOgulGT3KK7ZmZmZUb29vsOerHGNjYx00QbmX9tHA5RJJn7t0wT730husr0yPjJQQqCUBH9tRJ6GM240HoaylTxdaKcYrC8VJZhDwkkDthVKsYm+EPjs7q\/r6+ogovXTZ9haa8cr28uZqEKgigUYIpQleosiTJ0+2vmKMsopuWa0y8UquatmD0kCg3QQaJ5QasB1l1lUwfTRwu2+CNNdjfWUaSqSBQD0J+NiOOo1RpjGbnuzDOso0tJqbRk\/ukf1g3919S3NBUHMINIwAQllzg\/to4CqbhPdXVtk6lA0C5RDwsR0tPKIsB201cvXRwNUgF14KJvdU2TqUDQLlEPCxHUUoHXzBRwM7VK8jSZnc0xHsXBQCHSPgYzvaKKEMe9m0uT+t9pyoyUg+Grhjd4PDhZnc4wCLpBDwnICP7WhjhNIURL32U\/xNNlHYt2+fOnjwYLDDUNzho4F9uafYuccXS1FOCOQj4GM72gih1MtXtm7dqhYWFoJt9mSTBDnkf4cPH1bHjx9PfCWYjwbO59LtPVvv3CPvsGQmbHvZczUItIuAj+1oI4Ty6tWrgQ90dXWp4eHhFUIpRpONE6anpxHKdt0pMddhm7sKGIEiQKBEAghliXCLyDpsr1r7FWFx6z+1gaXrtqenJyiS\/l1E+chDKWbC4gUQqCeBxcXFoGLye3BwUJ07d86b9rMREaV2uzChlGhyaWlJTU5OBhGn\/dl0WS2U5nd79+5V8sNRHAE9E5Z3WBbHlJwg0GkCJ06cUPKjD4SygxYx95\/t7+9vCaAUKe7tJ7rIMmY5OjqqpqamVk3u0UJ57NgxtXHjxlZESVRZvMERy+KZkiMEOklAIkn5uXDhQiCYCGUnrRFz7bRCGTW5x8e+9YqaIlWxtFjuuvud6su77kx1DokgAIFqE\/CxHW1016teMrJt2zY1MDCg9Gd5IXXYC6h9NHC1b5nk0s298Ip6aO5FNbZ9kxrbflvyCaSAAAQqTcDHdrTRQineZG84YHfXmh7no4ErfcekLJwWS4kqJbrkgAAE\/CXgYzvaKKHM61o+GjhvnatyPmJZFUtQDgjkI+BjO4pQOtjcRwM7VK\/ySRHLypuIAkIgkYCP7ShCmWjWmwl8NLBD9bxIypilF2aikBCIJOBjO4pQOji0jwZ2qJ43SbVYMhvWG5NRUAi0CPjYjiKUDg7so4EdqudVUsTSK3NRWAgglE3xAYSyWpZmU4Jq2YPSQCANAR\/bUSLKNJZ9M42PBnaonpdJEUsvzUahG0zAx3YUoXRwWB8N7FA9b5PqjdTl1Vy8ostbM1LwhhDwsR1FKB2c00cDO1TP66QilrKDzy9fv66e3Hmnuv+ONV7Xh8JDoK4EfGxHEUoHb\/TRwA7Vq0VSeZ\/lcy+9EewNyy4+tTAplagZAR\/bUYTSwQl9NLBD9WqT9Oh3\/k4d\/c6VQCjZTL02ZqUiNSHgYzuKUDo4n48GdqherZLqST6MW9bKrFSmBgR8bEcRSgfH89HADtWrXVI9bildsTLJh3EQKDgPAAANfklEQVTL2pmYCnlIwMd2FKF0cDQfDexQvdomlUk+skHBfbevUfMP3VXbelIxCPhAwMd2FKF08CwfDexQvVonNbtimRVba1M3rnLScxJ1xP3PPkdmjGc53rX2FqfTrl5dVGMP\/Uf1g785q3p6epzO7VRihNKBvBbKb3\/72+q9731vcOZvf\/vb4Pett97K54rzMLti9UQf7If\/tuP+fe13f9hqabR4mcK0sHwt+L8pbObfK9NmEzSHpq70pP\/sH\/9B\/fQz9yKUpZPuwAW0UB46dEgNDw8HJfjGN74R\/P7Yxz7GZ094XHjtj9V\/urJGyUSff\/cnr6h71v0j9sN\/U9+\/y\/\/7D9UDD\/z7QNT+83\/9XnDexvf+WfD5F7\/4hZL\/v\/3t7wjW9P7mN78NPscd3f\/in9Rb33rjQfst\/+d\/Bb\/lfPFPyU+OD\/zbfx38\/tnP\/oda95Z\/Uvfee2\/w+Yc\/\/GHw+y8\/8sHg99\/+7bPB7z\/\/8w8Hv3\/y3\/8m+C3llePpp\/9L6OdPDe1MXX+zvXti5nSq\/M3rv\/76G+prX\/saEWWsV3j8Ty2U3\/zmN9XWrVuDmvzqV78Kfr\/jHe\/gs0c8fvdH\/1LJMhIZu5SG6r+N3h80TNizuf4sQif2F4ET\/5Ao74UXF4L7+nd\/9CfB91FdmeI7ckg3pP5bfr+7u6vV4ukuynVv+b\/q1lvf2lh\/u3Dhgvr0pz+tzp07R0TpsR5GFp0xyvpZVRq+HX\/9k6ABHNu+SY1tv61+laRGLYETO98QvGvBdzc\/r+7OFKGLEj79vRZFEKcn4GM7yhhlevsqHw3sUL1GJ9Wv7RIICKafrrBS+OKF0BbB++9YuyIiRADL8wEf21GE0sEffDSwQ\/Uan1Qa2rkX\/j7Y1UcaSpnwQ4RZLbcwxfAHP389KJx8J2tlzcMUQt0FShRYDVv62I4ilA6+46OBHapH0jcJmIKpI8xdd\/9pa+wJUOUS0OOAWvxEEG0xtMcEJSLUY4BsLFGuffLm7mM7ilA6WN1HAztUj6QWAVswpXGWCJPN1otxFTs6DBsvNCNDxLAY7p3Oxcd2FKF08BofDexQPZLGENAzZKUx192yRJnpXCZMEM2uUh0dys5JcmhBJDJMx9e3VD62owilg5f5aGCH6pE0BYGwKFMizKaLptldas4otQVRjxMSHaZwtpom8bEdbYRQXr58WQ0NDamlpaXA9UZGRtTY2FjLDc+cOaMmJiaCz0eOHFEDAwOhLuqjgWt6r1WiWlGiKWvn6to9m7W79MaEGretziphZApROAEf29HaC+Xy8nKwi44IY19fn9Kfd+7cGQiiiOjo6KiampoKHEL\/3dvbu8pBfDRw4V5uZLi4uKieeuop9dGPftSbhcNl8Tj\/s5+r\/zD+ZfX2bX\/ZWrOnu2jlmvfdvtaLt5dETaSROrh0l+IbKz0NHjd5+NiO1l4owxrGo0ePBl+LeEo0ef78eTU5Oam6urqU\/G\/Tpk2hUaWPBi5LGCRfeITf\/P\/vj\/9VICrSBSmbsdsCI2Nx5pKFdo3FaRHUk2ak9ElLLG6I\/MqxwzTRIb6x8s6DB0JZZltcSt6mUJp\/y8Xsz2YBcHZu\/iiHjPMNLVCyRvOGON0QTxEcc0s0PcNT0kSJkbklmqTTm2nf+Pvm7jL6b72ZdtjWa2GTaCQfGUfMK97cK9wrWe6VUhr8AjJtXESpxysfe+yxoCvWjiAlwrxy5cqKMUzNWd\/8s7Ozje9qFCbSnTQ4OKjgkY2FKWYinvq4EfGtfpuERKrmIW9gsI+b6wtv7jGqRTdsP9IC2pDQLPCNlVjgcZOHZsFer2XdfTnz1eOTIpB6Mo+LUIqBDxw4oGRTXw4IQAACEMhG4J577lFzc3PZTu7AWbWLKEX4Tp48GaDs7+9vjT2GiaSkcel61VGUCCYHBCAAAQhkIyAvbPblpc1Sw9oJZZjZ7JmuZhq7qzVuMk82l+AsCEAAAhDwmUDthfLatWtqfHxcbdiwIXTc0WV5iM+GpuwQgAAEIJCNQO2F0t5sQGMyu2XTbjiQDTFnQQACEICAzwRqL5Q+G4eyQwACEIBA5wkglJ23ASWAAAQgAIEKE0AoHYxjzqiVtYOyzKTuh9l1bXZXh9Xb7MKW\/2\/ZskVNT0+r7u7uumMK6qfHw3ft2tUI3zCNKmuMZbq\/3uEK\/7hBIGmf6TrfGGZ7KXNEZmZmVNjWoMKg6m0HQpnSU6UhEMNLw3\/p0qXW33UWAd3wb9u2Te3YsSOYFCV\/R20a3+QZw5rV\/Px8sAFDEx6i9K2jN+JIepBqmn8k7TOdsunxMpm9Nah8Pn36dOSDc9V9A6FM6YbmesumRA7mjGB5EoyLGprCJMxddNSwdetWtbCw0NqAP6VreZ1M7gt5OPjABz6gfv3rX0dGlE32D9PAcVtkeu0ICYW32xIzuQ++gVCm8E4zspJoyv6cIgsvk5hRtETO9mezUvbTs5cVzljoq1evBmfKpvrmm2oyZufVaeITEj3bEYRdiSb7B0J5ows66s1MPvgGQpmiWQp74ql6V0GKaiUmsSNIcfbDhw+r48ePrxp3DFuG07QuSB9u+ESjZ0yQJJT4x83xSr3PdEbUXp4m7aW8DzhsDNsH30AoU7gdQnnjFWRxQimiun\/\/\/taAvf05BWbvkyCUN19XZxuz6f4RtYWm906fogLyEPXEE09ETubxwTcQyhSGpuv1xszVuK5XG2NTuqfpfr5BICmibLJ\/IJLRIhnW\/Fax7UAoUwilJDG7Wn0YfE5ZrdhkdgSZZgmAzrApjBDKfEJZ96U0cftMF3GPVjmPpJmuUWWvYtuBUKb0NJaHxC8PcZ0OnhK7V8noeo3uem2ifyTtM+2VczsW1mXoxQffQCgdHIANB26+tkx3t5kvuTYXDSctMHbA7k1ShHKlUNpv5mmaf6TZZ9ob53YsqNlWmqfqCX6++QZC6egAJIcABCAAgWYRQCibZW9qCwEIQAACjgQQSkdgJIcABCAAgWYRQCibZW9qCwEIQAACjgQQSkdgJIcABCAAgWYRQCibZW9qCwEIQAACjgQQSkdgJIcABCAAgWYRQCibZW9qCwEIQAACjgQQSkdgJIcABCAAgWYRQCibZW9qCwEIQAACjgQQSkdgJIcABCAAgWYRQCibZW9qG0NAb2I9Pz8fmWpkZKT1v7GxsbbxlH1k9+3bpw4ePKh6e3vbdt0iLhT3dns7\/ya8EL0IpuTRXgIIZXt5czWPCLi8f7PsavksIC5C6fMDQdk+QP6dI4BQdo49V644gaoIpf1e0IpjW1U8F6GUk+03S\/hWX8pbPwIIZf1sSo0KIhAllBLdySFdr\/pdem9729uUvEJIDume3b17txoeHlYXL14MvtOvF5K\/7S7eLVu2qOnpadXd3R1acvN6OoH5yir57siRI2pgYKB1vvmao7BXnsWdb78eSuqju5m16En9Dh8+rJaWllbVzz7\/s5\/9rHrmmWfU1NRU0G0sXAcHB1tlNfOXL4kqC3JgsimMAEJZGEoyqhuBtEI5MTHREiotEuvXr2+Jn\/mmd2EkArpz586WsMW9CT7sHZd2ucyIraenR42PjysRRy1u9kt05XpPPPGEmpmZCYRLl\/mxxx4LTCgipoVdX7+vry\/IT6fdunWrmpycVF1dXcHDwunTp4P62vXT57\/66qvB9eQYHR1tiWZY\/fSDxLZt21aIf938i\/r4QwCh9MdWlLTNBNIKpRYJiQjDGnkzn0uXLqm5ubmWyJgRZpgwhHVbxglrVDetjkr37NkTCGmUCIVFr2b5X3vtNTU0NKREVEU85TDL+OMf\/1idP7\/yBc6mUMv5+\/fvb4l0lEnDytFm83M5CLQIIJQ4AwQiCKQVSlMYkoTy2WefVRKBhh1296ktQnq2q47CdLeueZ7drWlex+wSluhQC51OExXJmV2hdkRol\/Hs2bNBduaM4LCIV88strtddVl0l7aOWnFSCHSSAELZSfpcu9IEyhDKU6dOBeN6aQUgaSKMOdYogimHGeHagMO6OpOE0ixDXqE0l7aYom4LJkJZ6VujcYVDKBtnciqclkAZQikRZZyQ2WVLEko7AnvwwQfVoUOHIrs2k8b\/0nS9mmOMdkSZ1PUatgY0jDNdr2m9lHTtIIBQtoMy1\/CSQBlCKSBkMo+eHKOFxh73s6O8Xbt2tbpK7TFKU\/x27NgRjEHKoaNWHUXqCURx57\/nPe9JnMwTJ5Tr1q1bMVnJnswTNkZpi2KSmHvpTBTaawIIpdfmo\/BlEihDKGXCjz3GKHUwl4\/YdQpbV2gu\/5D0Ztdl2A5D9vhn3uUheqmHHVFKxGjXz14eYl+7v79\/RVd02ii6TNuTNwRMAggl\/gCBihPwfcMBV7xsOOBKjPRlE0AoyyZM\/hAogIDPW9i5VJ\/NBlxokbZdBBDKdpHmOhDIQaApAtKUB4IcrsCpHSCAUHYAOpeEAAQgAAF\/CCCU\/tiKkkIAAhCAQAcIIJQdgM4lIQABCEDAHwIIpT+2oqQQgAAEINABAghlB6BzSQhAAAIQ8IcAQumPrSgpBCAAAQh0gABC2QHoXBICEIAABPwhgFD6YytKCgEIQAACHSCAUHYAOpeEAAQgAAF\/CCCU\/tiKkkIAAhCAQAcI\/H\/gHlE+5vQEGAAAAABJRU5ErkJggg==","height":0,"width":0}}
%---
%[output:47d5323b]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"x_inf","rows":4,"type":"double","value":[["-0.0011"],["0"],["-0.0000"],["0"]]}}
%---
%[output:8ef96ea5]
%   data: {"dataType":"textualVariable","outputData":{"name":"u_inf","value":"0.3987"}}
%---
