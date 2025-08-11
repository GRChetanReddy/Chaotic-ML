%load('upo_states.mat')
N = 9;
a_vec = [ones(1, N^2), zeros(1, N^2)];

[xi_F, A, T] = sparrow_refine_UPO(xi_approx, T_approx, a_vec, N);

function [xi_F, A, T] = sparrow_refine_UPO(xi_approx, T_approx, a_vec, N)
    % Parameters from paper
    D_val = 1.3;
    alpha_val = -10;
    beta_val = 2;
    max_iter = 10;
    tol = 1e-4;
    
    % Project to Poincaré section
    b = 0;  % Paper uses b=0
    xi = xi_approx - a_vec'*(a_vec*xi_approx - b)/(a_vec*a_vec');
    T = T_approx;
    n = numel(xi);
    
    for iter = 1:max_iter
        % Integrate state + variational equations
        [xi1, A] = integrate_variational(xi, T, N, D_val, alpha_val, beta_val);
        
        % Newton system setup
        f_xi1 = network_dynamics(0, xi1, N);
        I_minus_A = eye(n) - A;
        M = [I_minus_A, -f_xi1; 
             a_vec,     0    ];
        rhs = [xi - xi1; 
               0];  % a·δξ=0
        
        % Regularized solve
        damping = 1e-3 * norm(M, 'fro');
        delta = (M + damping*eye(n+1)) \ rhs;
        
        % Apply updates
        xi = xi + delta(1:n);
        T = T + delta(end);
        
        % Check convergence
        if norm(delta) < tol
            break;
        end
    end
    xi_F = xi;
end

function J = analytical_jacobian(Z, N, D, alpha, beta)
    n = N^2;
    X = Z(1:n); Y = Z(n+1:end);
    J = sparse(2*n, 2*n);
    
    % 1. Add local nonlinear terms
    for j = 1:n
        idx = [j, j+n];  % Indices for oscillator j (real/imag)
        Z_j = [X(j); Y(j)];
        J(idx, idx) = local_jacobian(Z_j, D, alpha, beta);
    end
    
end

function J_local = local_jacobian(Z_j, D, alpha, beta)
    % Z_j = [X_j; Y_j]
    X = Z_j(1); Y = Z_j(2);
    
    % Jacobian blocks
    J11 = 1 - (3*X^2 + Y^2) - 2*beta*X*Y - 4*D;
    J12 = beta*(X^2 + 3*Y^2) - 2*X*Y + 4*alpha*D;
    J21 = -beta*(3*X^2 + Y^2) - 2*X*Y - 4*alpha*D;
    J22 = 1 - (X^2 + 3*Y^2) + 2*beta*X*Y -4*D;
    
    J_local = [J11, J12; J21, J22];
end

function [xi1, A] = integrate_variational(xi0, T, N, D, alpha, beta)
    n = numel(xi0);
    options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
    
    % Initial condition: [state; vec(identity matrix)]
    A0 = eye(n);
    Z0 = [xi0; A0(:)];
    
    % Integrate augmented system
    [~, Z] = ode45(@(t,Z) augmented_dynamics(t, Z, n, N, D, alpha, beta), [0 T], Z0, options);
    
    % Extract results
    Z_end = Z(end,:)';
    xi1 = Z_end(1:n);
    A = reshape(Z_end(n+1:end), n, n);
end

function dZ = augmented_dynamics(t, Z, n, N, D, alpha, beta)
    xi = Z(1:n);
    A = reshape(Z(n+1:end), [n, n]);
    
    J = analytical_jacobian(xi, N, D, alpha, beta);  
    dxi = network_dynamics(t, xi, N);
    dA = J * A;
    
    dZ = [dxi; dA(:)];
end

function dz = network_dynamics(t, Z, N)
    X = reshape(Z(1:N^2), N, N);
    Y = reshape(Z(N^2+1:end), N, N);
    D = 1.3;
    alpha = -10;
    beta = 2;
    
    Lx = zeros(N, N);
    Ly = zeros(N, N);
    % Interior points (2:N-1, 2:N-1)
    for i = 2:N-1
        for j = 2:N-1
            Lx(i,j) = X(i+1,j) + X(i-1,j) + X(i,j+1) + X(i,j-1) - 4*X(i,j);
            Ly(i,j) = Y(i+1,j) + Y(i-1,j) + Y(i,j+1) + Y(i,j-1) - 4*Y(i,j);
        end
    end
    
    % Zero-flux boundaries (Neumann)
    % Top/bottom edges (excluding corners)
    for j = 2:N-1
        Lx(1,j) = X(2,j) - X(1,j);     % Top edge
        Lx(N,j) = X(N-1,j) - X(N,j);   % Bottom edge
        Ly(1,j) = Y(2,j) - Y(1,j);
        Ly(N,j) = Y(N-1,j) - Y(N,j);
    end
    
    % Left/right edges (excluding corners)
    for i = 2:N-1
        Lx(i,1) = X(i,2) - X(i,1);     % Left edge
        Lx(i,N) = X(i,N-1) - X(i,N);   % Right edge
        Ly(i,1) = Y(i,2) - Y(i,1);
        Ly(i,N) = Y(i,N-1) - Y(i,N);
    end
    
    % Corners (average of adjacent edges)
    Lx(1,1) = (X(1,2) - X(1,1)) + (X(2,1) - X(1,1));  % Top-left
    Lx(1,N) = (X(1,N-1) - X(1,N)) + (X(2,N) - X(1,N)); % Top-right
    Lx(N,1) = (X(N,2) - X(N,1)) + (X(N-1,1) - X(N,1)); % Bottom-left
    Lx(N,N) = (X(N,N-1) - X(N,N)) + (X(N-1,N) - X(N,N));% Bottom-right
    
    dx = X - (X.^2 + Y.^2).*(X - beta*Y) + D*(Lx - alpha*Ly);
    dy = Y - (X.^2 + Y.^2).*(Y + beta*X) + D*(Ly + alpha*Lx);
    dz = [dx(:); dy(:)];
end