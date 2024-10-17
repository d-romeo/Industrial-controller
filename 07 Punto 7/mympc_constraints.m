% Constrained Model Predictive Control Function
% A: A matrix of the linear considered dynamic system
% B: B matrix of the linear considered dynamic system
% Q: status weight into the MPC cost function
% R: inputs weight into the MPC cost function
% S: final state weight (prediction horizon instant time) into the MPC cost function
% N: prediction horizon
% umin: inputs lower limit (scalar)
% umax: inputs upper limit (scalar)
% X: measured status at the current instant time
function u = mympc_constraints(A,B,Q,R,S,N,umin,umax,xmin,xmax,x,Ubar,Xbar,x0)
    m=size(B,2);
    n=size(A,1);
    
    %% Q and R matrices for Open-Loop MPC (with Q1)
    Qsig = blkdiag(kron(eye(N-1),Q),S);
    Rsig = kron(eye(N),R);

    %% A and B matrices for the Open-Loop MPC
    % A matrix
    Asig = A;
    for i = 2:N
        Asig = [Asig; A^i];
    end

    % B matrix
    Bsig = [];
    temp = [];
    for i = 1:N
        temp = zeros(size(B,1)*(i-1),size(B,2));
        for j = 0:N-i
            temp = [temp; A^(j)*B];
        end
        Bsig = [Bsig temp];
    end

    %% XBarsig, UBarsig
    Xbarsig = [repmat(Xbar,N, 1)]; 
    Ubarsig = [repmat(Ubar,N, 1)];
   
    %% H, F 
    H = Bsig'*Qsig*Bsig + Rsig;
    F = x'*Asig'*Qsig*Bsig - Xbarsig'*Qsig*Bsig - Ubarsig'*Rsig;
    f=F';

    %input and status constraints definition
    lb = [repmat(umin, N*m,1)];
    ub = [repmat(umax, N*m,1)];

    % AA e b per quadprog
    xmin = [repmat(xmin,N,1)];
    xmax = [repmat(xmax,N,1)];
    AA = [-Bsig;  
           Bsig]; 
    b = [-xmin+Asig*x; 
        xmax-Asig*x];

    options = optimset('Algorithm', 'interior-point-convex','Diagnostics','off', ...
        'Display','off');
    %solve the quadratic programming problem
    [U , FVAL, EXITFLAG]  = quadprog(H,f,AA,b,[],[],lb,ub,[],options);

    if (EXITFLAG ~= 1)
        U = quadprog(H,f,[],[],[],[],lb,ub,[],options);
    end
        

    %get the optimal input value (the receding horizon principle is applied)
    u = U(1:m);
end