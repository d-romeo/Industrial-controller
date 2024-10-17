clear all; close all; clc,

% A,B,C,D
A = [0 1 0 0 0 0 0 0 0 0 0 0;
    -20 -4 10 2 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0; 
    10 2 -20 -4 10 2 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0 0; 
    0 0 10 2 -20 -4 10 2 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0; 
    0 0 0 0 10 2 -20 -4 10 2 0 0;  
    0 0 0 0 0 0 0 0 0 1 0 0;  
    0 0 0 0 0 0 10 2 -20 -4 10 2;
    0 0 0 0 0 0 0 0 0 0 0 1; 
    0 0 0 0 0 0 0 0 10 2 -10 -2
    ];

B = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 1];

C = [1 0 0 0 0 0 0 0 0 0 0 0];

D = 0;

C_states = eye(12);
D_states = [0 0 0 0 0 0 0 0 0 0 0 0]';
%sampling time
Ts = 0.1;

%gain del sistema
sys = ss(A, B, C, D);
g = dcgain(sys);

%parametri di equilibrio
Ybar = 1;
Ubar = inv(g) * Ybar;
Xbar = -inv(A) * B * Ubar;

q = 10e-5;
r = 10e-6;

Q = [1+q 0 0 0 0 0 0 0 0 0 0 0;
    0 q 0 0 0 0 0 0 0 0 0 0;
    0 0 q 0 0 0 0 0 0 0 0 0;
    0 0 0 q 0 0 0 0 0 0 0 0;
    0 0 0 0 q 0 0 0 0 0 0 0;
    0 0 0 0 0 q 0 0 0 0 0 0;
    0 0 0 0 0 0 q 0 0 0 0 0;
    0 0 0 0 0 0 0 q 0 0 0 0;
    0 0 0 0 0 0 0 0 q 0 0 0;
    0 0 0 0 0 0 0 0 0 q 0 0;
    0 0 0 0 0 0 0 0 0 0 q 0;
    0 0 0 0 0 0 0 0 0 0 0 q;];
R = r;
%% 1.a) LQ CONTINUOUS TIME 
K = lqr(A, B, Q, R);
feedbackSystem = A - B*K;
poles_feedback = eig(feedbackSystem);
U_max = 15;
U_min = -15;
%open('LQ_continuous.slx')
%sim('LQ_continuous.slx')

%% 1.b) LQ DISCRETE TIME 
sysdis = c2d(sys,Ts); 
[A,B,C,D] = ssdata(sysdis);
g = dcgain(sysdis);
%pzmap(Gz); 
Ubar=inv(g)*Ybar;
Xbar_disc= (eye(12)- A)^-1*B*Ubar;
[K_disc, Pr,E]= dlqr(A, B, Q, R);
feedbackSystem = A - B*K_disc;
poles_feedback = eig(feedbackSystem);

%% 2) MPC controller
N = 10; 
x0= [0 0 0 0 0 0 0 0 0 0 0 0]'; 
S = Pr; 

%% 3) MPC controller: S = 0 e N variabile
N1 = 3;
N2 = 6;
N3 = 27; 


S =  [0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0;]; 
 
%% 4) MPC controller: control constraints
N1 = 3;
N2 = 6;
N3 = 27; 
umax = 15;
umin = -15;

%% 5) Constraints on the states
x0 = [9 0 9 0 9 0 9 0 9 0 9 0]'; 
xmax=[10 0 10 0 10 0 10 0 10 0 10 0]'; 
xmin= [-10 0 -10 0 -10 0 -10 0 -10 0 -10 0]'; 



%% 6) Kalman filter
x0 = [0 0 0 0 0 0 0 0 0 0 0 0]'; 
B_noise=[B ones(size(B,1),1)];
D_noise=[D_states zeros(size(D_states,1),1)];

power_vx=0.000001; 
power_vy=0.00001; 

% set Kalman filter parameters
QK=Q;               % variance of v_x
RK=10*QK(1);        % variance of v_y
NK=0;               % covariance v_x,v_y

x0K=1.5*x0;         %initial guess for x0

tsim=40;