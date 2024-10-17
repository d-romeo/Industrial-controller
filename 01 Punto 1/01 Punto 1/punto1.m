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
q1 = 10;
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
    0 0 0 0 0 0 0 0 0 0 0 q];

Q1 = [1+q1 0 0 0 0 0 0 0 0 0 0 0;
    0 q1 0 0 0 0 0 0 0 0 0 0;
    0 0 q1 0 0 0 0 0 0 0 0 0;
    0 0 0 q1 0 0 0 0 0 0 0 0;
    0 0 0 0 q1 0 0 0 0 0 0 0;
    0 0 0 0 0 q1 0 0 0 0 0 0;
    0 0 0 0 0 0 q1 0 0 0 0 0;
    0 0 0 0 0 0 0 q1 0 0 0 0;
    0 0 0 0 0 0 0 0 q1 0 0 0;
    0 0 0 0 0 0 0 0 0 q1 0 0;
    0 0 0 0 0 0 0 0 0 0 q1 0;
    0 0 0 0 0 0 0 0 0 0 0 q1];
R = r;
R1 = 10;
%% 1.a) LQ CONTINUOUS TIME 
K = lqr(A, B, Q, R);
feedbackSystem = A - B*K;
K1 = lqr(A, B, Q, R1);
feedbackSystem = A - B*K;
poles_feedback = eig(feedbackSystem);

N = 10; 
U_max = 15;
U_min = -15;


%% 1.b) DISCRETE TIME
sysdis = c2d(sys,Ts); 
[Ad,Bd,Cd,Dd] = ssdata(sysdis);
gd = dcgain(sysdis);
Ubar_disc=inv(gd)*Ybar;
Xbar_disc= (eye(12)- Ad)^-1*Bd*Ubar;
[K_disc, Pr,E]= dlqr(Ad, Bd, Q, R);
feedbackSystem = Ad - Bd*K_disc;
poles_feedback = eig(feedbackSystem);

open('sim_punto1.slx');
sim('sim_punto1.slx');

