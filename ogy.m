%load('refined_upo_state.mat')
N = 9;
f = network_dynamics(0, xi_F, N);
a_vec = [ones(1, N^2), zeros(1, N^2)];
n = length(xi_F);

denom = dot(f, f);
M = (eye(n) - (f*f')/denom)*A;
M_trans = M'; 

[eig_vecs, eig_vals] = eig(M_trans);

eig_vals = diag(eig_vals);
unstable_idx = find(abs(eig_vals) > 1); % Indices of unstable eigenvalues

u1_full = eig_vecs(:, unstable_idx(1));
u2_full = eig_vecs(:, unstable_idx(2));

if isreal(u1_full)
    u1 = real(u1_full);
    u2 = real(u2_full);
else
    u1 = real(u1_full);
    u2 = imag(u1_full);
end

thres = 0.1;
event = @(t, x) poincare_event(t, x, a_vec);

x0 = xi_F + 0.01*rand(size(xi_F));
x0 = x0 - a_vec' * (a_vec*x0)/(a_vec*a_vec');

t_total = 100*T;
t0 = 0;
sol.t = [];
sol.x = [];

while t0<t_total
    [t_seg, x_seg, te, xe, ie] = ode45(@(t, x) network_dynamics(t, x, N), ...
        [t0, t_total], x0, odeset('Events', event, 'RelTol', 1e-5, 'AbsTol', 1e-7));

    sol.t = [sol.t; t_seg];
    sol.x = [sol.x; x_seg];

    if ~isempty(te)
        xe = xe(end, :)';
        dist = norm(xe - xi_F);

        if dist < thres
            A_mat = [dot(u1, u1), dot(u1, u2); 
                     dot(u2, u1), dot(u2, u2)];
            rhs = [dot(u1, xi_F - xe);
                   dot(u2, xi_F - xe)];
            epsilon = A_mat \ rhs;

            % Apply perturbation and project to section
            x_perturbed = xe + epsilon(1)*u1 + epsilon(2)*u2;
            x_perturbed = x_perturbed - a_vec' * (a_vec * x_perturbed) / (a_vec * a_vec');
            x0 = x_perturbed;

        else
            x0 = xe;
        end
        t0 = te;

    else
        break;
    end
end

plot(sol.t, sol.x(:, (5-1)*N + 5)); % Re(W_{5, 5})
xlabel('Time');
ylabel('Re(W_{5,5})');
title('Controlled UPO');



%% Vizualise the phase space of refined upo %%
% tspan = [0 5000];
% Z = xi_F;
% [t, Z] = ode45(@(t, Z) network_dynamics(t, Z, N), tspan, Z);
% 
% x = Z(:, 1:N^2);
% y = Z(:, N^2+1:end);
% 
% x = x(:, 41);
% y = y(:, 41);
% 
% figure;
% plot(x, y);
% xlabel('x');
% ylabel('y');
% title('Phase space of oscillator (5,5)');

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

function [value, isterminal, direction] = poincare_event(t, x, a_vec)
    value = a_vec * x; % a_vec · ξ = 0
    direction = 1;     % Trigger when crossing from negative to positive
    isterminal = 1;    % Stop integration at event
end