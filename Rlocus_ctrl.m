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

A = [0 1 0 0;... %[output:group:14413a11] %[output:28ac8dce]
    A1(1,1) 0 A1(1,2) 0;... %[output:28ac8dce]
    0 0 0 1;... %[output:28ac8dce]
    A1(2,1) 0 A1(2,2) 0] %[output:group:14413a11] %[output:28ac8dce]
B = -[0 ; B1(1); 0 ; B1(2)] %[output:8c0fde09]
C = [0 0 1 0] %[output:2a3a9b65]
D = [0] %[output:84889453]
%%
[num,den] = ss2tf(A,B,C,D) %[output:60988226] %[output:06e878d2]
%%
sys = tf(num,den);
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
s = tf('s');

% your plant
G = tf(num, den);

% a placeholder PID (tune next)
Kp = 300;
Ki = 10;
Kd = 5;

C = pid(Kp, Ki, Kd);

% closed-loop transfer (unity feedback)
T = feedback(C*G, 1);
[y_imp, t_imp] = impulse(T);
y_imp = y_imp/4;

info = impulseinfo(t_imp, y_imp) %[output:49f89039]

plot(t_imp, y_imp) %[output:199741e0]
grid on %[output:199741e0]
title('Impulse Response') %[output:199741e0]


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
%[output:28ac8dce]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A","rows":4,"type":"double","value":[["0","1.0000","0","0"],["-11.8839","0","-2.5891","0"],["0","0","0","1.0000"],["47.6759","0","38.6973","0"]]}}
%---
%[output:8c0fde09]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"B","rows":4,"type":"double","value":[["0"],["-0.8402"],["0"],["4.2977"]]}}
%---
%[output:2a3a9b65]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"C","rows":1,"type":"double","value":[["0","0","1","0"]]}}
%---
%[output:84889453]
%   data: {"dataType":"textualVariable","outputData":{"name":"D","value":"0"}}
%---
%[output:60988226]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"num","rows":1,"type":"double","value":[["0","0","4.2977","0","11.0147"]]}}
%---
%[output:06e878d2]
%   data: {"dataType":"matrix","outputData":{"columns":5,"name":"den","rows":1,"type":"double","value":[["1.0000","-0.0000","-26.8134","-0.0000","-336.4381"]]}}
%---
%[output:49f89039]
%   data: {"dataType":"textualVariable","outputData":{"header":"struct with fields:","name":"info","value":"            Peak: 7.3151\n    SettlingTime: 0.5410\n       Overshoot: 279.9462\n"}}
%---
%[output:199741e0]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQvYVVW19wehCQUqiCWKHDlgSZ+JxAFvaEejQ6KYmRcwHyjNK5rXR0s9Shwz8E6A4h1QuXiDxOzgBaygjqiBx0JNRI9RmPIAat+RSuV7xvJbr5uXvd+x9h77subcv\/U8PMD7zrn3Wr\/5H2P+95xjrd1u06ZNm4QDAhCAAAQgAAEIBESgHQYmoNHiVCEAAQhAAAIQSAhgYBACBCAAAQhAAALBEcDABDdknDAEIAABCEAAAhgYNAABCEAAAhCAQHAEMDDBDRknDAEIQAACEIAABgYNQAACEIAABCAQHAEMTHBDxglDAAIQgAAEIICBQQMQgMBmBE444QT59a9\/LXfeead8+ctfhg4EIACBXBLAwORyWDipZiNwzDHHyLPPPisTJ06U4cOHN\/Tya21gFi9eLKNGjdriGrt27Spf+MIX5Mwzz5RBgwY1lAFvDgEI5J8ABib\/Y8QZNgGBZjQwn\/70p+Vf\/\/Vfk9HVB4KvWrVKXnzxRdl6661l\/vz58rnPfa4JRp5LhAAEKiWAgamUHP0gUEUCrQ3MscceK88880yyjXPjjTfK8uXLZY899pApU6bIz3\/+c5k6daq0b99exowZI6NHj07O5LjjjpOnn35arr32Wpk9e3bSZ8cdd5RLLrlEDj300KTNV77yFXn11VeT1\/j85z+f\/Ozf\/u3fZOXKlS0\/a70C8+abb8qPf\/xjWbJkibzzzjuy8847y7e+9S056aSTWgg89NBDyXnqa2+77bbyzW9+Uy644ALZaquttqCUrsD06tVLnnjiic1+n77397\/\/fTnllFPk\/fffT1alHnzwQXnrrbdk1113TVZovvGNbyT99PfXXXed6PuvXbtWtttuO\/nqV78qP\/jBD0QNkrJQZsrpT3\/6kzz++ONJP2V2\/vnnt7z30qVLk7b\/\/d\/\/nfxszz33lLPPPlsGDx6c\/D99nUsvvTR5nfvuu0+22WYbOf744+Xcc8\/NdC7aqBxOVZQXLwWBKAlgYKIcVi4qNAKtDYxOjP\/1X\/8lvXv3Fv2dTnwrVqxITMxnPvMZ2XvvvWXSpEnJZf7yl7+UXXbZJZlMtY\/+Xg3AunXrEuPxiU98QtQ0qJmpxMB85zvfkSeffDIxDj169JC5c+eKTvjXXHONHHXUUcl76nvruWmbp556Su666y4577zzkv+3PkoZGF2FUROmxm3s2LHJNtMNN9wgP\/nJT2TYsGFy2GGHJa+r76cGYsCAAXLLLbfI+PHjk98ddNBBsmzZssS8qYG6+uqrE0bXX399YjbUxPzzP\/9z8tpqdrTvkCFDkj76vsrpxBNPTP6+9dZb5cMPP5Sf\/exnyUpQ+jpq3vbbbz9R86Xn9fe\/\/13uvvtu2X\/\/\/c1zKZdTaBrmfCFQbwIYmHoT5\/0gUIRAawOTrkSkJmDhwoXy3e9+Vz75yU8mqyydO3eWkSNHJmZBV2N0FSXtoysfZ5xxRvIuJ598crLKkRqCSgyMrkLoKoyahn79+snbb7+drLT80z\/9k3Tp0iU5Lz0\/PQ8t+tWJP\/1bzUgpA6MrJAceeGDy6w8++CBZBdLXVZPwyCOPJNe41157yV\/\/+lf57W9\/K5\/61KeSbSY1K2poJk+enKyiqKHS61MTpatS+p76GmrqtI2u0Og1zJgxI3mv22+\/XX70ox+1vEbKKF310TZpvyOPPDLpn\/5fTZNy0ENXtmbNmpWswJx11lnmuZTLiUCBAATaJoCBQSEQyAGBUgYmXSV44YUXkom7cNtFzc28efPkqquukqOPPrrFwNx2221yyCGHJFc1YcIEufnmmxOTcfHFF1e0AqMrLbo9pMcOO+wg++yzjxx++OHyta99LfnZwQcfLP\/zP\/9TlKKubui2TuFRqohX2+hWl5oRXS1as2aNHHDAAUVfV1d71OQsWrQouTZdvVGD86UvfSm5xhEjRiSrLqnxSK9fX+wXv\/iF6KqSmiPlp++h7zVnzhwZOHBg8n7pOWpR8cMPP9zyOt\/+9rflsssuS9rcdNNNySpP+trWuZTLKQey5BQgkGsCGJhcDw8n1ywEShmY9Fbml156KZnc+\/TpI48++miCRVdatDaktYFJTY+2+Y\/\/+I+kjubUU0+Viy66qMXA6NZI3759NzMgaV1MsbuQdELX+hG9U+r3v\/990u\/yyy9PaknULL322mty5ZVXttTVpOOmBkCNRDEDU2jGZs6cKVpfoj9TY6J93njjjWRrRuto1FwUHlroq3UqevzhD39Iin515UUNk27r6IqUrgilBkbrdXTFRA+9Dq2v0dUkXb3R99D3KjQwv\/rVr5Jr0y28xx57rOV1Uo76OrrNpFt0heaorXMpl1OzaJ\/rhEClBDAwlZKjHwSqSKCaBkZrPdIC1bSw94orrki2WLT49bnnnmvZdtLCWF2B0GLYYgZGJ\/eXX35Z2rVr12J4tL++jv5O6z\/SLRg1MLryoashWpejqzWpybAMjG47HXHEEUmdjxbP6h99HTUZuoWk22BqbtavXy+\/+93vkloc\/b9uOf3lL3+RfffdN3kLLTLWIt4NGzaIrlrpypFuARVu\/WhNjNa06CqS1rGk56+Fv\/pvPbRwWP+03kJqy8BY56J99TqycqqivHgpCERJAAMT5bByUaERqKaBUeOgZkXNiRa06taK3kGkWznpioxO\/rrSoqZFt3+0bboqU7gC8y\/\/8i9JcazWqGidh9a8aAHvvffe27LykBan6vuefvrpiQnRlaG0TqX1WJQq4tXVHeWgqysLFiyQ3XbbrcVIfPGLX0wKc\/UctQZo3Lhxyfnrlo5uCZ122mmy++67Jyspajx0dUm3h9IVmE6dOsnQoUOlZ8+eyV1Jukozffr0pAZH31eNl76vrszo7+64447EQOlr6Gulr9OWgbHOpVxOoWmY84VAvQlgYOpNnPeDQBEC1TQwutqidy3pbdTdu3dPaja0\/kIPXcHQrSR90q6aES1cVTOiWyYPPPCA9O\/fv6WWJt2+0hUPrfXQlZeNGzcm9Sm6eqHFq1pUrIdO9GoM1AypUdJ6HX2fjh07bnG1bd1Gndb1qLFQg5HeRq2Fs3ruWpirWzvpreP6My3I1RUfLS7WW7i1YFe31wqLeLXmRVdn1LB16NAhuQVczVZ6KA9dmdHbqLUQWGtpdBVLeeiRxcBY51IuJwIFAhBomwAGBoVAIBICtX6CboiYihmPEK+Dc4YABLYkgIFBFRCIhAAGZsuBxMBEIm4uAwJFCDSlgfnHP\/6RFNLp3rkWJ+qSuC6z6x0THBAIlQAGBgMTqnY5bwhUQqApDYx+KlPzog+20j16\/VuflaGFjhwQgAAEIAABCOSfQFMaGH28ud6NcOGFFyYjpLc\/6sOvtOhRiwA5IAABCEAAAhDIN4GmNDD6wKpp06Yl36vSrVu35A4DfT6DPtSKAwIQgAAEIACB\/BNoSgOjz3fQZ1rosyb0+1j0tkp97kOpGpj0O2fyP5ycIQQgAAEIQGBzAvr1H\/q9XbEdTWlg9CFYf\/zjH+Xaa69Ntoz02RD6PTG6CtO1a9ctxli\/wVa\/RI6jOAH4tK0M+MDHmzvQEBryaChW\/TSlgdEHXeljw\/VhW+mhX+ymjxfXb9FtfcQ6+J6AKOwLH5KrR0vox6YHI2LMVknpFrHqpykNjD4uXL8sTldg9Emi+hTOE088MVmB0ad3YmDKC5VYg6M8Cs2XPOBTLQL26xBjGBhbJc2Xg5rSwOj3vug2kt51pN90q3Uw55xzjgwZMqSoAkgebYfO2LFjRf9wFCcAH\/TjjQ00hIY8Gop1DmtKA1OuEGId\/HI5lGqvt6HrlwNyFCcAn7aVAR87cmCEhmyVsALjYRRtXwwMycMjbiYf9OPRj\/ZFQ2jIo6FY5zBWYDKoItbBz3DpmZqQXEmumYRSohH6senBiBizVcIKjIdRtH0xMCQPj7iZfNCPRz+swNj0iLG2GcU6h7ECY8eGxDr4GS49UxOSBxN0JqGwAlMxJmKMGKtYPCLRzmEYmAyqwMCQPDLIpGQTJh\/049EPKzA2PWKMFRhbJU3aAgPDBOSRPskV\/Xj0g4Gx6RFjGBhbJU3aAgPDBOSRPskV\/Xj0g4Gx6RFjGBhbJU3aAgPDBOSRPskV\/Xj0g4Gx6RFjGBhbJU3aAgPDBOSRPskV\/Xj0g4Gx6RFjGBhbJU3aohwD8\/q6jbLklQ2yeOV60X\/37NpBLhraK\/k71oPkwQTt0Tb6senBiBizVVK6RTlzmOd96t2Xu5AyEC9n8CcseDUxLoP7dJFdu3SQJa+sl1lPvyHLL90vwzuF2YTkSnL1KBf92PRgRIzZKsHAeBhF2zergVGjogamtVkZM+uFxNTMH9M\/SkYkV5KrR9jox6YHI2LMVgkGxsMo2r5ZDIwalCNuXCaTR\/SVwX2234LF8CnLkm2kKSP7RseJ5Epy9Yga\/dj0YESM2SrBwHgYRdvXMjBqXnSVRY2L1ruUOtTE6O+LGZyQ4ZFcSa4e\/aIfmx6MiDFbJRgYD6No+1oGZvHKDTLr6TXm6opuMc1cuia6rSSSK8nVE\/zox6YHI2LMVgkGxsMo2r6WgdGVleMHdZeRA3cyGex9xW9KbjOZnXPagORKcvVIE\/3Y9GBEjNkqwcB4GEXbty0Dk9a+ZL3LSFdh9BbrmGphSK4kV0\/wox+bHoyIMVslGBgPo2j7tmVgKtkW0lWYh87oH82zYUiuJFdP8KMfmx6MiDFbJRgYD6No+7ZlYCopzE2fFRPLKgzJleTqCX70Y9ODETFmqwQD42EUbd9SBqbc7aNCQDGtwpBcSa6e4Ec\/Nj0YEWO2SjAwHkbR9i1lYDz1LB\/ddt0lU+Fv3sGSXEmuHo2iH5sejIgxWyUYGA+jaPuWMjCVbB+lkDzmJ2+gSa4kV48m0Y9ND0bEmK0SDIyHUbR9ixkYz\/ZRCqrreYuSrx0I\/YseSa4kV0\/wox+bHoyIMVslGBgPo2j7FjMw1SjEjWUbieRKcvUEP\/qx6cGIGLNVgoHxMIq2bzED49k+im0bieRKcvUEP\/qx6cGIGLNVgoHxMIq2bzEDU60n6uo20rrrDg6aHcmV5OoRMPqx6cGIGLNVgoHxMIq2b2sDU436lxRWDNtIJFeSqyf40Y9ND0bEmK0SDIyHUbR9WxsY\/fJGrYGZP6a\/+5pjuBuJ5Epy9QQC+rHpwYgYs1WCgfEwirZvawOj5kWPi4b2qso1h76NRHIluXoCAf3Y9GBEjNkqwcB4GEXbt7WBqUYBbyGscr7NOo+QSa4kV48u0Y9ND0bEmK0SDIyHUbR9Cw1MNetfUmChbyORXEmunuBHPzY9GBFjtkowMB5G0fYtNDDVrH8pNDAzl66pSk1NIwaB5Epy9egO\/dj0YESM2SrBwHgYRdu30MBUu\/4lhRbyU3lJriRXT\/CjH5sejIgxWyUYGA+jaPsWGphq17+k0EK+nZrkSnL1BD\/6senBiBizVYKB8TCKtm9qYGpR\/1K4jbR45XqZMrJvcBxJriRXj2jRj00PRsSYrRIMjIdRtH1TA1OL+pcUmr72ETcuC\/KpvCRXkqsn+NGPTQ9GxJitEgyMh1FQfRcuXChXXnml\/OUvf5Fdd91VLrvsMtl3332LXkNqYGpV\/5K+qX49wUNn9A\/u26lJriRXT\/CjH5sejIgxWyUYGA+jYPquWrVKjjzySJk6dWpiWu677z6ZN2+ezJw5U9q1a7fFdaQGptZ1KrV+\/VoNEMmV5OrRFvqx6cGIGLNVgoHxMAqmr668vPPOOzJ+\/PhM55wamFoV8KYnEerzYEiuJNdMgVSiEfqx6cGIGLNVgoHxMAqm74gRI2SvvfaSl156SV555RXZZZdd5OKLL5Z+\/fqV3ELSX2w48nbZft5JMmrUKBk9enRNrnfApNfk2bN2q8lr1+pFV69eLT169KjVywf\/uvBpewjhY0scRmjIVsnHLaZPny4zZszYrIvuPMR2tNu0adOm2C7Kup5hw4bJP\/7xD7nzzjule\/fucscdd8gtt9wiixYtkk6dOm3RXVdgnnxmRVJku\/zS\/ayXd\/0+xDoYPh3y6dAjevRj04MRMWarhBUYD6Ng+o4cOVL22WcfOeecc1rOWVdfbrrpJtl\/\/\/2LGpgZj\/62at9A3RaoEOtgSK4kV0\/wox+bHoyIMVslGBgPo2D6XnLJJdK+fXsZN27cZgZGV2HU2LQ+dAXm1JueSH5crW+gLgUrxDoYkivJ1RP86MemByNizFYJBsbDKJi+zz\/\/vJxwwgly9913y5577inTpk2Tm2++OdlC6tixY1ED83\/OfyAxL4P7bF\/T69TnwZw5+4Wab1VV8yJIriRXj57Qj00PRsSYrRIMjIdRUH3nzJkjkyZNknfffVf69OkjP\/zhDxMzU+zQFZh6GRh9\/9C+F4nkSnL1BD\/6senBiBizVYKB8TCKtu9uew6S7UdMrNuqiN6uffyg7jJy4E5BMCW5klw9QkU\/Nj0YEWO2SjAwHkbR9q23gdFC3p5dO9S83qZaA0ZyJbl6tIR+bHowIsZslWBgPIyi7dtz0KHSb\/SVMn9M\/7pcY2iFvCRXkqsnMNCPTQ9GxJitEgyMh1G0fXcedq6cffb36roionUw6647OAimJFeSq0eo6MemByNizFYJBsbDKNq+nxk1VR68bETN70AqBBjSA+1IriRXT\/CjH5sejIgxWyUYGA+jaPs2wsCE9EA7kivJ1RP86MemByNizFYJBsbDKMq+r6\/bKF+6eL6snXpMXa8vpDoYkivJ1RMc6MemByNizFYJBsbDKMq++mC5I695rO4GRt93woJX61Y47Bk8kivJFf14CNh9iTFizFYJBsbDKMq+aiSOGjdb3pxxWt2vL5RCXpIrydUTHOjHpgcjYsxWCQbGwyjKvrqFdNBBB8lrv1ta9+vTQl79+oK8P9CO5Epy9QQH+rHpwYgYs1WCgfEwiravfpXAqlWr6n59oTzQjuRKcvUEB\/qx6cGIGLNVgoHxMIq2b6MMTCiFvCRXkqsn+NGPTQ9GxJitEgyMh1G0fRtlYBRoCHUwJFeSqyf40Y9ND0bEmK0SDIyHUbR9G2lgQnigHcmV5OoJfvRj04MRMWarBAPjYRRt30YamBAeaEdyJbl6gh\/92PRgRIzZKsHAeBhF27fRBibv30xNciW5eoIf\/dj0YESM2SrBwHgYRdu3kQYmhEJekivJ1RP86MemByNizFYJBsbDKNq+jTQwCjXvhbwkV5KrJ\/jRj00PRsSYrRIMjIdRtH0bbWDyXshLciW5eoIf\/dj0YESM2SrBwHgYRdu30QYm74W8JFeSqyf40Y9ND0bEmK0SDIyHUbR9G21g9Esd9dCvFcjjQXIluXp0iX5sejAixmyVYGA8jKLt22gDk\/dCXpIrydUT\/OjHpgcjYsxWCQbGwyjavo02MAo2z4W8JFeSqyf40Y9ND0bEmK0SDIyHUbR982Bg8lzIS3IluXqCH\/3Y9GBEjNkqwcB4GEXbNw8GJs+FvCRXkqsn+NGPTQ9GxJitEgyMh1G0ffNgYPJcyEtyJbl6gh\/92PRgRIzZKsHAeBhF2zcPBibPhbwkV5KrJ\/jRj00PRsSYrRIMjIdRtH3zYGAUbl4LeUmuJFdP8KMfmx6MiDFbJRgYD6No++bFwOS1kJfkSnL1BD\/6senBiBizVYKB8TCKtm+eDIw+zG7kwJ1yxZrkSnL1CBL92PRgRIzZKsHAeBhF2zcvBiavdyKRXEmunuBHPzY9GBFjtkowMB5G0fbNi4HJayEvyZXk6gl+9GPTgxExZqsEA+NhFG3fvBiYxSs3yJmzX5Dll+6XK9YkV5KrR5Dox6YHI2LMVgkGxsMo2r55MTAKOI93IpFcSa6e4Ec\/Nj0YEWO2SjAwHkbR9s2TgcnjnUgkV5KrJ\/jRj00PRsSYrRIMjIdRtH3zZGDyWMhLciW5eoIf\/dj0YESM2SrBwHgYRds3TwYmj4W8JFeSqyf40Y9ND0bEmK0SDIyHUZB958+fL2effbbMmTNHBg4cWPQaMDAkD4+4mXzQj0c\/2hcNoSGPhvI0h3muo3Xfdps2bdpUzRcM6bXeeOMNOe644+Ttt9+WW2+9NQgDo3zzVshLciW5euIe\/dj0YESM2SphBcbDKLi+o0ePluHDh8sNN9wg119\/fTAGJm+FvCRXkqsn+NGPTQ9GxJitEgyMh1FQfWfMmCGLFy+WW265RQYPHmwamMKLGzVqlKj5adQx9vG1MmCXDjK8b6dGncJm77t69Wrp0aNHLs4ljycBn7ZHBT62amGEhmyVfNxi+vTponNc4bFq1apyXiKItk25haSfZtSEzJ07V7p165bJwORp8CcseDURl34vUh4OPh3y6dCjQ\/Rj04MRMWarhBUYD6Mg+n7wwQdy7LHHyumnny5DhgxJzjnLCkyeDEze7kQiuZJcPcGPfmx6MCLGbJVgYDyMguj70ksvJQamc+fOLee7Zs0a6dq1q5x88slyyimnbHEdeazgzlMhL8mV5OoJfvRj04MRMWarBAPjYRRs39BWYBR0ngp5Sa4kV0\/wox+bHoyIMVslGBgPo2D7hmpgtAZm5MCdGs6d5Epy9YgQ\/dj0YESM2SrBwHgYRds3j1tIefpKAZIrydUT\/OjHpgcjYsxWCQbGwyjavnk0MHkq5CW5klw9wY9+bHowIsZslWBgPIyi7ZtHA7N45QbR26nnj+nfcO4kV5KrR4Tox6YHI2LMVgkGxsMo2r55NDAKOy93IpFcSa6e4Ec\/Nj0YEWO2SjAwHkbR9s2rgcnLnUgkV5KrJ\/jRj00PRsSYrRIMjIdRtH3zamDyUshLciW5eoIf\/dj0YESM2SrBwHgYRds3rwZGa2BeX7dRpozs21D2JFeSq0eA6MemByNizFYJBsbDKNq+eTUwebkTieRKcvUEP\/qx6cGIGLNVgoHxMIq2b14NjALPQyEvyZXk6gl+9GPTgxExZqsEA+NhFG3fPBuYPBTyklxJrp7gRz82PRgRY7ZKMDAeRtH2zbOBGT5lmRw\/qHtDv1KA5Epy9QQ\/+rHpwYgYs1WCgfEwirZvng1MHu5EIrmSXD3Bj35sejAixmyVYGA8jKLtm2cDk4dCXpIrydUT\/OjHpgcjYsxWCQbGwyjavnk2MHn4SgGSK8nVE\/zox6YHI2LMVgkGxsMo2r55NjAKvdF3IpFcSa6e4Ec\/Nj0YEWO2SjAwHkbR9s27gWn0nUgkV5KrJ\/jRj00PRsSYrRIMjIdRtH3zbmAaXchLciW5eoIf\/dj0YESM2SrBwHgYRds37wam0V8pQHIluXqCH\/3Y9GBEjNkqwcB4GEXbN+8GptF3IpFcSa6e4Ec\/Nj0YEWO2SjAwHkbR9s27gVHwjSzkJbmSXD3Bj35sejAixmyVYGA8jKLtG4KBaWQhL8mV5OoJfvRj04MRMWarBAPjYRRt3xAMTCO\/UoDkSnL1BD\/6senBiBizVYKB8TCKtm8IBqaRdyKRXEmunuBHPzY9GBFjtkowMB5G0fYNwcA0spCX5Epy9QQ\/+rHpwYgYs1WCgfEwirZvCAamkV8pQHIluXqCH\/3Y9GBEjNkqwcB4GEXbNwQDo\/AbdScSyZXk6gl+9GPTgxExZqsEA+NhFG3fUAxMo+5EIrmSXD3Bj35sejAixmyVYGA8jKLtG4qBaVQhL8mV5OoJfvRj04MRMWarBAPjYRRt31AMTKO+UoDkSnL1BD\/6senBiBizVYKB8TCKtm8oBqZRdyKRXEmunuBHPzY9GBFjtkowMB5G0fYNxcDoADSikJfkSnL1BD\/6senBiBizVYKB8TCKtm9IBqYRhbwkV5KrJ\/jRj00PRsSYrRIMjIdRtH1DMjCN+EoBkivJ1RP86MemByNizFYJBsbDKNq+IRmYRtyJRHIluXqCH\/3Y9GBEjNkqwcB4GEXbNyQD04hCXpIrydUT\/OjHpgcjYsxWCQbGwyjaviEZmEZ8pQDJleTqCX70Y9ODETFmqwQD42EUTN\/3339frrrqKvnpT38qGzdulN69e8vll18u\/fr1K3oNIRkYvYB634lEciW5eoIf\/dj0YESM2SrBwHgYBdN38uTJ8sgjj8i0adOkW7duMnHiRLn\/\/vtlyZIlURiYet+JRHIluXqCH\/3Y9GBEjNkqwcB4GAXTd\/HixdK5c+eWFRdNDocccoisWLFCOnbsuMV1hLYCU+9CXpIrydUT\/OjHpgcjYsxWCQbGwyjYvlOnTpWFCxfKvffeG8UKTL0LeUmuJFdP8KMfmx6MiDFbJRgYD6Mg+86bN0\/Gjx8vs2bNkl69epU0MIW\/GDVqlIwePTq31zv\/hb\/Ks3\/aKGOHdKvLOa5evVp69OhRl\/cK8U3g0\/aowcdWNYzQkK2Sj1tMnz5dZsyYsVmXVatWlfMSQbRtt2nTpk1BnGmVT\/LDDz+Ua665RhYtWiQ333yz9OzZs+Q7hLaFpBdSz0JePh3y6dATnujHpgcjYsxWCSswHkbB9FXPduGFF8ratWtl0qRJ0qlTpzbPPUQDU89CXpIrydUT\/OjHpgcjYsxWCQbGwyiYvjNnzhT9M3fuXNl6663N8w7RwNSzkJfkSnI1g6iNBujHpgcjYsxWCQbGwyiYvsOGDZOXX35Z2rdvv9k533PPPTJgwIAtriNUA9Ozawe5aGjxup5qDhbJleTq0RP6senBiBizVYKB8TCKtm+IBqaedyKRXEmunuBHPzY9GBFjtkowMB5G0fYN0cC8vm6jHHHjMll+6X41HxeSK8nVIzL0Y9ODETFmqwQD42EUbd8QDYwORr3uRCK5klw9wY9+bHowIsZslWBgPIyi7RuqganXnUgkV5KrJ\/jRj00PRsSYrRIMjIdRtH1DNTD1uhOJ5Epy9QQ\/+rHpwYgYs1WCgfEwirZvqAamXoW8JFeSqyf40Y9ND0bEmK0SDIyHUbR9QzUwi1dukAkLXpX5Y\/rXdGxIriRXj8DQj00PRsSYrRIMjIdRtH1DNTA6IPUo5CXqMOJ7AAAgAElEQVS5klw9wY9+bHowIsZslWBgPIyi7RuygalHIS\/JleTqCX70Y9ODETFmqwQD42EUbd+QDUw9CnlJriRXT\/CjH5sejIgxWyUYGA+jaPuGbGDqUchLciW5eoIf\/dj0YESM2SrBwHgYRdsXA0Py8IibyQf9ePSjfdEQGvJoKOQ5rK3rbrdp06ZNHjDN0Df0wa91IS\/JleTqyQPox6YHI2LMVgkrMB5G0fYN3cDUupCX5Epy9QQ\/+rHpwYgYs1WCgfEwirZv6Aam1oW8JFeSqyf40Y9ND0bEmK0SDIyHUbR9Qzcw+jA7PS4a2qsmY0RyJbl6hIV+bHowIsZslWBgPIyi7Ru6gan1nUgkV5KrJ\/jRj00PRsSYrRIMjIdRtH1DNzA6MLUs5CW5Zkuur6\/bKPonPXp27SD6p9kP9GMrAEbZYswm2ZwtYpjDio0cdyFl0HMMg1\/LQl6Sa3ERqVmZ9fQaefx3b8hb733UZtcuHxuWP67\/yNCoiRk5cKfk97Xa5ssg84Y1QT82ehhhYGyVsALjYRRt3xgMTC0LeUmuH0tfDcmSVzbIzKVrRA2KGpPen94oRw\/uWzI+UqPz0d9vJIZGjUxqaqINrP9\/YejHHmEYYWBslWBgPIyi7RuDgallIS\/J9SPpq\/lQzgf03l5GDuwug\/tsn\/y8XD76GvpaeqiJiX1Vplw+0SaaNi4MRhgYj+5jmMOKXT9bSBlUEcPg17KQt9mTa7qCoownj+jbYlxSaVXKR19XzYyu6MRsZCrlkyF0o2kCIwyMR8wxzGEYmAoVEMvg16qQt5mTq5qMj7bnti+5UuLlk76HbkkVM0gVyjo33bx8cnMhNTwRGGFgPPKKZQ5rzYAVmAyqiGXwa1XI26zJNd3qsepVqsVn8coNcsSNy5LVmCkjS9fUZJB0rppUi0+uLqrKJwMjDIxHUrHMYRiYClQQy+DXqpC3GZNrWu\/y0Bn9zVuhq80nNU6xrMZUm08FIZ77LjDCwHhEGsschoGpQAWxDH6tCnmbLbmmBiKLeVG51YKPbivpaozelj1\/TP8KVJ2fLrXgk5+rq86ZwAgD41FSLHMYBqYCFcQy+LUq5G2m5Jqal+WX7pdZSbXkE8NqTC35ZB6knDeEEQbGI9FY5jAMTAUqiGnwa1HI2yzJtdyVl1RqteaT1sZcNHS3IG+5rjWfCkI+d11ghIHxiDKmOayQA0W8GVQR0+DXopC3GZJrOTUvrSVVDz6Fdypl3drKIP26NKkHn7pcSA3fBEYYGI+8YprDMDBlKiGmwa9FIW\/syTWtN6m0aLaefNRo6RiriUkfpFem3OvevJ586n5xVXpDGGFgPFKKaQ7DwJSphJgGX7dBdEKu5m24MSfX1Lx4HiRXbz7pOesTgas5zmWGTebm9eaT+cRy1BBGGBiPHGOawzAwZSohpsGvRSFvzMlVVzP08BiBRvEZPmVZ8n1Med9SahSfMtNAQ5vDCAPjEWBMcxgGpkwlxDb41S7kjTW5pnUv5dxxVExajeQTwpZSI\/mUmQoa1hxGGBiP+GKbw1IWFPFmUEVsg1\/tQt4Yk2u6DVON1YtG88n7llKj+WRIAQ1vAiMMjEeEsc1hGJgy1BDb4Fe7kDe25Jre0XP8oO7JY\/u9R1746LjrF0NWWozs5VCqf1741Or6qvG6MMLAeHQU2xyGgSlDDbENfrXrYGJLrtUudM4Tn3RLKU\/PjMkTnzLSQl2bwggD4xFcbHMYBqYMNcQ2+DqJzVy6pmqPoI8puVZz6yiVWN74pCtMen55+BqCvPEpIzXUrSmMMDAescU2hzW1gdm0aZNcd911MmfOHNm4caP069dPfvzjH0uPHj2KaiTGwa9mIW8sybXaW0d5NTDpeeXlawhi0U\/r5KF60kP\/1rvBPvr3ey0\/K2yfti38Wc+uHVr+27ndRunSZfuW\/\/fs2jH5Hiw9tF1hW89EF2rfWDVUrfGIcQ5TNk1ZxHvPPffIbbfdJnfddZd89rOflauvvlqWLl0q8+bNaxoDU81C3liSR7VXpvJuYPT89GsIzpz9gjTymTEh60eNR2pQFq9cn\/xbD601Sk2FGo303x\/\/3XGzXJOakcIfpqZHf\/bcK2s2MzDp++rvtF36vq3fs\/D90vOI0eyErKFqmZS2XgcDUw\/KdXqPY445Rr7+9a\/LCSeckLzje++9J3vvvbc8\/PDDsvvuu29xFjEOfjULeWNIHt6n7bYl3RD4NHI1JgQ+6UqKGpPUqKQmJTUGg\/t0aVkVqfZTkLMyar3qU7jikxot\/Ts2o5OVT52mmNy9TYxzWNOuwKhZmTp1quy7774tQvvKV74i5513nhx22GFNYWCqWcgbQ\/LQh77ppHPR0F5VTz6h8ElXY3RCrmdtTB75FK6saL1YalZ0pSo1KtU2KfU0wZUYndbbVa23sdLfVz2AMrxgORpqvV1X+P\/CVa\/U\/KXmtfA00j6aL+qpgwwommYXoWkNzOc\/\/3mZPXu29O\/fv2Wwhw0bJieeeKIcffTRRQ1M4Q9HjRolo0ePrlRLuek3YNJr8uxZu7nPZ\/Xq1SXrh9wvXocXePZPG2Xs42tl\/ujiNVDeUwiNz81PbZCHX\/yrHL5HJzl1n4\/rLrwcSvXPA58\/v\/O+rHn3fXlm9UZRPei\/9RiwS4fkz\/C+nWp1+Zlet1GMlIseKQ\/9f+HPiv1e2++87VYt19W988f\/Lvx52qDw9xaM9Dxat3v1rf8r22yzzWY\/Lmybnmfrc9P\/lzq\/ts47vQ5tU+yarOuo9e+nT58uM2bM2OxtVq1aVeu3rfvrN2UNjK7ATJo0SQ488MAW4IMHD5ZLLrlEDj300KIGJsbBr1YdTDmffuqucOMN08LdWn6SCpFP4bdb1\/q5MY3io9eoKyu6wpJ+8tYVlpEDu+fuU3WjGFUar4WrGqVWONLXLlzpsN5PV32KHVv\/bb3stFP3zX5VWOsTY92Pxarw92whlUMr521HjBghumV08sknJ2e6fv16GTRokDz66KPSq9eWWwixDn616mBCS66F8qxV4W7he4TMJ91W0uuplZGpF5\/UsGgNi467TmrpllA1HlhYy7RXL0a1vIZavjZ82qYb6xzWlCsw999\/v0ycODG5C6l79+4ybtw4WblyZXJbdbEj1sGvVh1MqMmjloW7sRiY9DrSIl\/9v65WVXPCr5V+0k\/+s55ek9xtpassjaph8U7etWLkPa+89IcPBiYvWqzLeaiBufvuu5M7kAYOHJg8B2annYo\/Nh4D0\/aQhJo8qvFN01nEGiqfYtdWaGTUxFSj6LmafIqtsuh5HtC7S+62hbJoJ21TTUblvG8obeGDgQlFq3U\/z1gNjIKsxgPtQkwetXjibilhhsjHCjLlpysbExa8lqxqeL43ysMn5G0hi3Fsq3jlXG+5bT0aKve9Qmwf6xzWlFtI5Qow1sFXDtUo5A0xeejqi9ZAVGMFwdJTiHysa0p\/nxoZ3aLRgtiPtmi2L2vFIyuf1Kxo0Wf6foV1LHr7dwi3tGZli4HJTiqrhrK\/YlwtY53DMDAZdBrr4OulV6OQN7Tkka6+LL90vwyj728SGp9Krzh9dsqSV9ZvZjDSepn0uSGtn0qrfNpv9\/EdJB8\/cO295AmzqXHRfoUPjatmHU6l11yvfs2ioUp5woctpEq1E32\/mA1MNQp5Q0se+tA6z5ZHuYIPjU+511eqfWERrbYpfDhc+v+0b+FtroWPvE9NT6wrK1lZN6uG4JOVAAamOqQifJWYDYwOl7cOJqTkmt4WXK\/VF+UbEp9GhC98bOowapsRfDAwdhQ1aYvYDYy3Diak5KGrL7V8aF2xEAmJTyNCHD42dRhhYGyVlG4R6xxGDUwGVcQ6+Omle+tgQkmujVh9YQXGDrBQ9GNfSe1awAgD41FXrHMYBiaDKmId\/PTSvXUwoSTXRqy+YGDsAAtFP\/aV1K4FjDAwHnXFOodhYDKoItbBL7x0Tx1MCMm1Hl8ZUEpKIfDJEAY1awIfGy2MMDC2SthC8jCKtm8zGBhPHUzek2s9vrCxLfHnnU+jAxc+9gjACANjqwQD42EUbd9mMDCeOpi8J9dGrr6whWSnhbzrx76C2reAEQbGo7JY5zC2kDKoItbBL7x0Tx1MnpNro1dfMDB2gOVZP\/bZ16cFjDAwHqXFOodhYDKoItbBb33pldbB5Dm5Nnr1BQNjB1ie9WOffX1awAgD41FarHMYBiaDKmId\/NaXXmkdTF6Tax5WXzAwdoDlVT\/2mdevBYwwMB61xTqHYWAyqCLWwW996ZXWweQ1ueZh9QUDYwdYXvVjn3n9WsAIA+NRW6xzGAYmgypiHfzWl15pHUwek2teVl8wMHaA5VE\/9lnXtwWMMDAexcU6h2FgMqgi1sEvdumV1MHkMbnqU3cnLHhV5o\/pn2GEa9skj3xqe8XlvTp8bF4wwsDYKindItY5DAOTQRWxDn6xS6+kDiaPybVRT90txjSPfDLIvm5N4GOjhhEGxlYJBsbDKNq+zWRgKqmDyVtybdR3HpUKgLzxyVugwsceERhhYGyVYGA8jKLt20wGppI6mLwl1zytvmhQ5I1P3gIVPvaIwAgDY6sEA+NhFG3fZjIwOohaB7P80v2kZ9cOmcY0T8k1b6svGBhbQnnSj322jWkBIwyMR3mxzmHUwGRQRayDX+rSy91GylNyzdvqCwbGDrA86cc+28a0gBEGxqO8WOcwDEwGVcQ6+KUuvdxtpLwk1zyuvmBg7ADLi37sM21cCxhhYDzqi3UOw8BkUEWsg9\/WpZdzO3VekmseV18wMHaA5UU\/9pk2rgWMMDAe9cU6h2FgMqgi1sFv69LLuZ06D8k1r6svGBg7wPKgH\/ssG9sCRhgYjwJjncMwMBlUEevgt3Xp5dTB5CG56urL8YO6y8iBO2UY0fo2yQOf+l5xee8GH5sXjDAwtkpKt4h1DsPAZFBFrIPf1qWXUwfT6OSa59UXVmDsAGu0fuwzbHwLGGFgPCqMdQ7DwGRQRayDb1161tupG51c87z6goGxVMZzcmxCMLIYNToHWefX6N\/HOodhYDIoK9bBty496zZSI5NH3ldfMDCWypicbUIwshg1MgdZ55aH38c6h2FgMqgr1sG3Lj3rNlIjk0feV18wMJbKmJxtQjCyGDUyB1nnloffxzqHYWAyqCvWwc9w6Zmeytuo5BHC6gsGxlZZo\/Rjn1l+WsCIGhiPGmOdwzAwGVQR6+BnuHTJso3UqOSa1+e+tObaKD5ZxjcPbeBjjwKMMDC2Skq3iHUOw8BkUEWsg5\/h0iXLNlIjkmsoqy+swNgqa4R+7LPKVwsYYWA8iox1DsPAZFBFrIOf4dKTJtZTeRuRXENZfcHA2CprhH7ss8pXCxhhYDyKjHUOw8BkUEWsg5\/h0pMm1jZSvZOrrgpNWPBq8o3ZIRz15hMCk8JzhI89YjDCwNgqYQvJwyjavs1uYKxtpHom19fXbUwM1UVDe8ngPtsHobl68gkCSKuThI89ajDCwNgqwcB4GEXbt9kNTLqNpCsePbt22GKc65lc1UzNXLpG5o\/pH4ze6sknGCgFJwofe9RghIGxVYKB8TAKpu\/7778vV111lfz0pz+VjRs3Su\/eveXyyy+Xfv36Fb0GDMxH20hqXnTlo\/VRr+Sqqy9H3LhMJo\/oG8zqi7KqF59gApAVmLKHCg1hYMoWTUGHWOewpqyBmTx5sjzyyCMybdo06datm0ycOFHuv\/9+WbJkCQamRJToXT9qHtZdd3DDDIy1leUJ8Fr2ZfJh8vHqCw2hIY+GMDAeejnru3jxYuncuXPLiosmh0MOOURWrFghHTt23OJsYx38coel1J0\/9Uiu6erLQ2f0L7qNVe611LN9PfjU83qq\/V7wsYnCCANjq6R0i1jnsKZcgWk9zFOnTpWFCxfKvffeywpMG1FSagWkHslVt7D0mDKyryeOG9K3HnwacmFVelP42CBhhIGxVYKB8TAKsu+8efNk\/PjxMmvWLOnVa8v6Dr0oda+Fx6hRo2T06NFBXq\/3pAdMek3mj+4hO2+7VctLrV69Wnr06OF96ZL9\/\/zO+3Lq3DeS9w3xqDWfEJkUnjN87BGEUduM4LM5n+nTp8uMGTM2++GqVatsoQXWIvoVGP3kctRRR7UMy6OPPio77rijfPjhh3LNNdfIokWL5Oabb5aePXuWHLpYl98q0WqxYt5afzoM4Qsb22JZaz6VjGOe+sDHHg0YsQJjq4QVGA+jXPb94IMPZN26dS3ntsMOO0i7du3kwgsvlLVr18qkSZOkU6dObZ47BuZjPMWKeWuZXEN7aF0xIdWSTy6DrsyTgo8NDEYYGFslGBgPo2D6zpw5U\/TP3LlzZeuttzbPGwOzOSJdEdGHyKW3VNcquYZ623RrQdWKjyncQBrAxx4oGGFgbJVgYDyMguk7bNgwefnll6V9+\/abnfM999wjAwYM2OI6MDCbI1FjsfcVv2m5pbpWybW1UQpGYK1OtFZ8QuWBwSt\/5NAQBqZ81XzcI9Y5LPoaGM+gp31jHXwPm8JamFok1xi2jlK+teDjGbu89YWPPSIwwsDYKmEFxsMo2r4YmC2HtnAVptrJNcTvO2pL\/NXmE1ugwcceURhhYGyVYGA8jKLti4EpPrTps1ku2LdDyVvQyxVFal4Ka2zKfY28tWfyYfLxahINoSGPhmKdw9hCyqCKWAc\/w6WbTbqet0huOWonOXpwdR4wN2HBq6J3OoX0ZY0WJCYfJh9LI9bv0RAasjTS1u9jncMwMBlUEevgZ7h0s4majdPufl5+N\/ZAs63VIK17CfHrAthCska39O+ZnG12MMLA2CphC8nDKNq+GJi2h\/ar1\/5GPrfz9q7H\/Mdyy3QxUkw+TD7e5IiG0JBHQ7HOYazAZFBFrIOf4dIzNdHk+r1HNiTPhdHalXKPGOteChkw+TD5lBsTrdujITTk0VCscxgGJoMqYh38DJeeqYkm1z990EXOnP2CTB7RtywTk668jBy4U8uD8TK9aUCNmHyYfLxyRUNoyKOhWOcwDEwGVcQ6+BkuPVOTNLmmXzOg3xithsQ6tL2anpjNizJg8mHysWLB+j0aQkOWRtr6faxzGAYmgypiHfwMl56pSWFyzbKiom30bqMlr2yQ2Ap2iwFj8mHyyRRIbTRCQ2jIo6FY5zAMTAZVxDr4GS49U5PWyTU1KHpX0QG9t0+2lHp27Sivr3tP9HdqXGJfdSkEx+TD5JMpkDAwFWMixtpGF+schoHJEDKxDn6GS8\/UpFTyULMy6+k1iWlJj8F9uiSmpmfXDpleO4ZGJFcMjFfHaAgNeTQU6xyGgcmgilgHP8OlZ2pCciW5ZhJKiUbox6YHI2LMVknpFrHOYRiYDKqIdfAzXHqmJiRXkmsmoWBgKsZEjBFjFYtHRGKdwzAwGVQR6+BnuPRMTUiuJNdMQsHAVIyJGCPGKhYPBsaDLvy+GBiSh0fFTD7ox6Mf7YuG0JBHQ7HOYazAZFBFrIOf4dIzNSG5klwzCYUVmIoxEWPEWMXiYQXGgy78vhgYkodHxUw+6MejH1ZgbHrEWNuMYp3DWIGxYyPaAqgMl56pCcmDCTqTUFiBqRgTMUaMVSweVmA86MLvG6t7rdbIkFxJrh4toR+bHoyIMVslpVvEOoexApNBFbEOfoZLz9SE5EpyzSQUVmAqxkSMEWMVi4cVGA+68PtiYEgeHhUz+aAfj360LxpCQx4NxTqHsQKTQRWxDn6GS8\/UhORKcs0kFFZgKsZEjBFjFYuHFRgPuvD7YmBIHh4VM\/mgH49+WIGx6RFjbTOKdQ5jBcaODe5CMhiRPJigM4RRySbox6YHI2LMVknpFhgYD73A+8Y6+NUaFpIrydWjJfRj04MRMWarBAPjYRRtXwwMycMjbiYf9OPRD1tINj1ijC0kWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gIDwwTkkT7JFf149IOBsekRYxgYWyVN2gID0\/bAjx07VvQPR3EC8EE\/3thAQ2jIo6FY57Cm\/zLH+fPny9lnny1z5syRgQMHFtVIrIPvCYjCvvBpzk8\/6KdaBOzXIcaIMVslpVvEqp+mNjBvvPGGHHfccfL222\/LrbfeioGpMEJiDY4KcWzRDT5MPl4toSE05NFQrPppagMzevRoGT58uNxwww1y\/fXXY2AqjJBYg6NCHBiYMsGhHxsYjDAwtkpYgfEwCqrvjBkzZPHixXLLLbfI4MGD2zQwI0eOlKeeeiqo6+NkIQABCEAAAkpgn332kVmzZkUHoylXYLRifdSoUTJ37lzp1q2baWCiG3UuCAIQgAAEIBA4gegNjJqVo446qmWY\/vM\/\/1POOOMMOf3002XIkCHJz60VmMDHmNOHAAQgAAEIREcgegPzwQcfyLp161oGbu3atTJixAjp3Llzy8\/WrFkjXbt2lZNPPllOOeWU6AaZC4IABCAAAQjERiB6A5NlwFiByUKJNhCAAAQgAIH8EMDAsIWUHzVyJhCAAAQgAIGMBDAwGUHRDAIQgAAEIACB\/BDAwORnLDgTCEAAAhCAAAQyEsDAZARFMwhAAAIQgAAE8kMAA1NiLNavXy\/\/\/u\/\/Lr\/85S\/lE5\/4hHzta1+TH\/7wh7LNNtvkZ\/TqdCY\/+9nP5Nprr5W33npLunfvLpdddlly63mx4+9\/\/7v85Cc\/kZtuuknuuOMO+fKXv1yns2zs25TD6P7770\/46B1x+hyiM888U77xjW809gJq\/O7l8NGHTN5+++2iMfjZz35Wvve97yVPzI75KIdPyuG1116TYcOGJY+FUA3FfmRl9MILL8hhhx0mn\/zkJzdD8vTTT29292lsvLLy0evWx4v84Ac\/kBUrVkinTp3k1FNPFX0yfWgHBqbEiJ111lnyt7\/9Ta677jrZtGlTcov1XnvtJRdffHFoY+w63z\/84Q9y5JFHJk8sPuCAA+Sxxx6TCy64QJ544gnZcccdt3jtr3\/963LQQQfJzJkzE3bNYGDKYaRPdP7ud78rOkn3799ffvWrX8lJJ50kCxYskF69ernGKq+dy+Hz+OOPJzE2e\/Zs0cfn6\/\/1mU2qu9122y2vl+g6r3L4pG+kj4c49thjk0dEfPOb34zewJTDaOnSpXLuuefKkiVLXOMSUudy+OiHzKFDh8rxxx8v3\/nOd+TFF1+U8847T+68807ZZZddQrpswcAUGa733ntP+vXrJ+pod99996TFr3\/96+Rbq9XFN9OhKy\/q1idPntxy2foFmPqJ+IQTTtgCxfLly2XvvfeWQYMGydVXX90UBqYcRppoXn\/99ZaHKCpAfcz3+PHj5eCDD45SWuXw0U+E77zzjuy7774tLNTo6feVxWqGy+GTQrnxxhtl1apVyX\/V2MW+AlMOo0cffTT5apif\/\/znUcZTsYsqh49+KPjRj34kixYtCp4PBqbIEOoS5OGHHy4vv\/xysn2kx5tvvpkk1WeeeSZ56F2zHKeddprssccecs4557Rcsn5C3nrrrZMttVJHMxmYShkpu2XLliVLt7\/4xS+kS5cuUcqqUj76QeLBBx+UKVOmJCtUhQ+fjAlUuXw0P40ZM0bmzZsn48aNawoDUw6j++67TyZOnJisaL700kvJSrHy0u22WI9y+OiHgd\/\/\/veyww47JN8HuO222yY7DCFuY2Ngiij62WefTZbXVPzp8e677yarMloT06NHj1jjYIvr0u+M2m+\/\/ZJl\/PTQpKmfkq+55hoMjEjyvVqVMNKJSLeTLr30Ujn00EOj1VQlfCZMmCBTp05NlrR1K1INcaxHOXx0W1snGtXM\/vvvn2znNsMKTDmMdFLWVZhvfetbyQr6woULkxWqe+65RwYMGBCljMrho\/lbt2hvu+22REO6ra0fovTLHnW1M6QDA1NiBUaLwHRvMC0EW716dVLboeYm1k\/KxYSrzl5rES688MKWX2vS1MKvsWPHYmBEpBJGurx9xRVXJFtHBx54YEg5o+xzrYSPvonu1Wty1dU\/3TLRrbYYj3L4XHnllQmCtBavWQxMOYyKaUS\/Ikbz2Pe\/\/\/0YJVRWDtLVqd\/85jeJiUkPLeLVlXatHQrpwMAUGS39lKMFuw888IDsueeeSYu0uFALxJrpSJcbb7311pbLVnOnn250larU0UxbSOUy0iVuvQtJ77SJtXC3UBfl8NFPz3qn38CBA1teQpe3+\/TpIxdddFGUoVcOH737T28qaNeuXcJCi3i32mqrxATrVlusRzmMtGZPGalhSQ8tVtWcHtoEnXU8y+Ezf\/785E5RLYxPDzV4X\/ziF0VvXgnpwMCUGC39ZLNhw4Zk+VoNjSZRTRLnn39+SOPrPldNBlqwq8lRk+dDDz0kl19+eUvNxnPPPZck0dYFqM1kYMphpEW8WgStHHfddVf3+ITwAuXwUaN81113JXdp6daIFvWqUdaiQzXOMR7l8Gl9\/c2yAlMOI71jUreL7r777mS7X++Y1BqYuXPnyhe+8IUYJZTcaJE1T+t8prlctaN3sumNF\/qB9N577235wB4KJAxMiZHSGg\/dZ37yySeTQl69lVj\/r592mu3Q\/WTd6vjzn\/+cfKrR4t30E7LWKjz\/\/PNJstDtNQ0EPXT5X1kpO13+LSwCjpFfVka6\/zxt2rQtnlGhd7gV1hnFxigrnw8\/\/DD50KAFqlp3tt122yV3u8X+LfFZ+TSrgdHrzspINaT1eVoA\/r\/\/+7\/JBwVdeRkyZEhsYbXZ9WTlo5305gF9npfeEanPotL8c5TuVO4AAAErSURBVMQRRwTHBwMT3JBxwhCAAAQgAAEIYGDQAAQgAAEIQAACwRHAwAQ3ZJwwBCAAAQhAAAIYGDQAAQhAAAIQgEBwBDAwwQ0ZJwwBCEAAAhCAAAYGDUAAAhCAAAQgEBwBDExwQ8YJQwACEIAABCCAgUEDEIAABCAAAQgERwADE9yQccIQgAAEIAABCGBg0AAEIAABCEAAAsERwMAEN2ScMAQgAAEIQAACGBg0AAEIQAACEIBAcAQwMMENGScMAQhAAAIQgAAGBg1AAAIQgAAEIBAcAQxMcEPGCUMAAhCAAAQggIFBAxCAAAQgAAEIBEcAAxPckHHCEIAABCAAAQhgYNAABCAAAQhAAALBEcDABDdknDAEIAABCEAAAhgYNAABCEAAAhCAQHAE\/h92fOnPp5XTtAAAAABJRU5ErkJggg==","height":337,"width":560}}
%---
