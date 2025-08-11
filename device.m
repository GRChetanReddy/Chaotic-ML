N = 9;
load('C3_upo.mat');

%% Input of pattern plus %%
% I = zeros(N, N);
% for i = 1:9
%     I(5, i) = 1;
%     I(i, 5) = 1;
% end

%% Input pattern for cross %%
I = eye(N);
% + fliplr(eye(N));
% I(5, 5) = 1;


Z = [xi_F; xi_F; zeros(2*N^2, 1); reshape(I, N^2, 1)];

f = network_dynamics(0, Z, N);
a_vec = [ones(1, N^2), zeros(1, N^2)];
n = length(xi_F);
f = f(1:2*N^2);

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
event = @(t, x) poincare_event(t, x, a_vec, N);

x0 = xi_F + 0.01*rand(size(xi_F));
x0 = x0 - a_vec' * (a_vec*x0)/(a_vec*a_vec');

Z0 = [x0; x0; zeros(2*N^2, 1); reshape(I, N^2, 1)];

sol.t = [];
sol.z = [];
sol.w = [];

t_total = 100*T;
t0 = 0;

while t0<t_total
    [t_seg, Z_seg, te, Ze, ie] = ode45(@(t, Z) network_dynamics(t, Z, N), ...
        [t0, t_total], Z0, odeset('Events', event, 'RelTol', 1e-5, 'AbsTol', 1e-7));

    sol.t = [t_seg; sol.t];
    sol.z = [Z_seg(1:2*N^2); sol.z];
    sol.w = [Z_seg(2*N^2+1:4*N^2); sol.w];

    if ~isempty(te)
        xe = Ze(1, 2*N^2)';
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
            Z0 = x_perturbed;

        else
            Z0 = Ze;
        end
        t0 = te;

    else
        break;
    end
end

function dz = network_dynamics(t, Z, N)
    Z_X = reshape(Z(1:N^2), N, N);
    Z_Y = reshape(Z(N^2+1:2*N^2), N, N);
    W_X = reshape(Z(2*N^2+1:3*N^2), N, N);
    W_Y = reshape(Z(3*N^2+1:4*N^2), N, N);
    Z_Px = reshape(Z(4*N^2+1:5*N^2), N, N);
    Z_Py = reshape(Z(5*N^2+1:6*N^2), N, N);
    I = reshape(Z(6*N^2+1:end), N, N);

    D = 1.3;
    alpha = -10;
    beta = 2;
    gamma = 2.44;
    
    Lzx = zeros(N, N);
    Lzy = zeros(N, N);
    Lwx = zeros(N, N);
    Lwy = zeros(N, N);
    % Interior points (2:N-1, 2:N-1)
    for i = 2:N-1
        for j = 2:N-1
            Lzx(i,j) = Z_X(i+1,j) + Z_X(i-1,j) + Z_X(i,j+1) + Z_X(i,j-1) - 4*Z_X(i,j);
            Lzy(i,j) = Z_Y(i+1,j) + Z_Y(i-1,j) + Z_Y(i,j+1) + Z_Y(i,j-1) - 4*Z_Y(i,j);
            Lwx(i,j) = W_X(i+1,j) + W_X(i-1,j) + W_X(i,j+1) + W_X(i,j-1) - 4*W_X(i,j);
            Lwy(i,j) = W_Y(i+1,j) + W_Y(i-1,j) + W_Y(i,j+1) + W_Y(i,j-1) - 4*W_Y(i,j);
        end
    end
    
    % Zero-flux boundaries (Neumann)
    % Top/bottom edges (excluding corners)
    for j = 2:N-1
        Lzx(1,j) = Z_X(2,j) - Z_X(1,j);     % Top edge
        Lzx(N,j) = Z_X(N-1,j) - Z_X(N,j);   % Bottom edge
        Lzy(1,j) = Z_Y(2,j) - Z_Y(1,j);
        Lzy(N,j) = Z_Y(N-1,j) - Z_Y(N,j);
        Lwx(1,j) = W_X(2,j) - W_X(1,j);     
        Lwx(N,j) = W_X(N-1,j) - W_X(N,j);   
        Lwy(1,j) = W_Y(2,j) - W_Y(1,j);
        Lwy(N,j) = W_Y(N-1,j) - W_Y(N,j);
    end
    
    % Left/right edges (excluding corners)
    for i = 2:N-1
        Lzx(i,1) = Z_X(i,2) - Z_X(i,1);     % Left edge
        Lzx(i,N) = Z_X(i,N-1) - Z_X(i,N);   % Right edge
        Lzy(i,1) = Z_Y(i,2) - Z_Y(i,1);
        Lzy(i,N) = Z_Y(i,N-1) - Z_Y(i,N);
        Lwx(i,1) = W_X(i,2) - W_X(i,1);     
        Lwx(i,N) = W_X(i,N-1) - W_X(i,N);   
        Lwy(i,1) = W_Y(i,2) - W_Y(i,1);
        Lwy(i,N) = W_Y(i,N-1) - W_Y(i,N);
    end
    
    % Corners (average of adjacent edges)
    Lzx(1,1) = (Z_X(1,2) - Z_X(1,1)) + (Z_X(2,1) - Z_X(1,1));  % Top-left
    Lzx(1,N) = (Z_X(1,N-1) - Z_X(1,N)) + (Z_X(2,N) - Z_X(1,N)); % Top-right
    Lzx(N,1) = (Z_X(N,2) - Z_X(N,1)) + (Z_X(N-1,1) - Z_X(N,1)); % Bottom-left
    Lzx(N,N) = (Z_X(N,N-1) - Z_X(N,N)) + (Z_X(N-1,N) - Z_X(N,N));% Bottom-right
    Lwx(1,1) = (W_X(1,2) - W_X(1,1)) + (W_X(2,1) - W_X(1,1));  % Top-left
    Lwx(1,N) = (W_X(1,N-1) - W_X(1,N)) + (W_X(2,N) - W_X(1,N)); % Top-right
    Lwx(N,1) = (W_X(N,2) - W_X(N,1)) + (W_X(N-1,1) - W_X(N,1)); % Bottom-left
    Lwx(N,N) = (W_X(N,N-1) - W_X(N,N)) + (W_X(N-1,N) - W_X(N,N));% Bottom-right
    
    dzx = Z_X - (Z_X.^2 + Z_Y.^2).*(Z_X - beta*Z_Y) + D*(Lzx - alpha*Lzy) + Z_Px;
    dzy = Z_Y - (Z_X.^2 + Z_Y.^2).*(Z_Y + beta*Z_X) + D*(Lzy + alpha*Lzx) + Z_Py;

    if any(Z_Px, "all") || any(Z_Py, "all")
        dwx = W_X - (W_X.^2 + W_Y.^2).*(W_X - beta*W_Y) + D*(Lwx - alpha*Lwy) + gamma*I.*(Z_X - W_X);
        dwy = W_Y - (W_X.^2 + W_Y.^2).*(W_Y + beta*W_X) + D*(Lwy + alpha*Lwx) + gamma*I.*(Z_Y - W_Y);
    else
        dwx = W_X - (W_X.^2 + W_Y.^2).*(W_X - beta*W_Y) + D*(Lwx - alpha*Lwy);
        dwy = W_Y - (W_X.^2 + W_Y.^2).*(W_Y + beta*W_X) + D*(Lwy + alpha*Lwx);
    end

    dz = [dzx(:); dzy(:); dwx(:); dwy(:); Z_Px(:);  Z_Py(:);  I(:)];
end

function [value, isterminal, direction] = poincare_event(t, x, a_vec, N)
    x = x(1:2*N^2);
    value = a_vec * x; % a_vec · ξ = 0
    direction = 1;     % Trigger when crossing from negative to positive
    isterminal = 1;    % Stop integration at event
end
