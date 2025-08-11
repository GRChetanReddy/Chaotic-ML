N = 9;
tspan = [0 5000];

%% Random Initialization %%
% psi = pi/4;
% X = 1*cos(psi)*ones(N, N) + 1*randn(N, N);
% Y = 1*sin(psi)*ones(N, N) + 1*randn(N, N);
% Z = [X(:); Y(:)];

%% Initialization for standing waves%%
Z = initialize_system(1, 1, 9);

%% Required for power spectrum %%
[t, Z] = ode45(@(t, Z) network_dynamics(t, Z, N), tspan, Z);

x = Z(:, 1:N^2);
y = Z(:, N^2+1:end);

x = x(:, 41);
y = y(:, 41);

[crossing_states, crossing_times] = record_poincare_crossings(N, tspan);
[xi_approx, T_approx] = find_upo_candidate(crossing_states, crossing_times);


%% vizualize the poincare section at (5,5) %%
x = crossing_states(41, :);
y = crossing_states(N^2+1+41, :);

figure;
scatter(x, y);
xlabel('x');
ylabel('y');
title('crossing points of the poincare section oscillator (5,5)');

%% Visualize the power spectrum %%
% transient = t < 900;
% x = x(~transient);
% t = t(~transient);
% t_uniform = t(1):dt:t(end);
% x = interp1(t, x, t_uniform, 'pchip');
% 
% n = length(x);
% window = hann(n);
% x_windowed = x .* window;
% y1 = fft(x_windowed);
% P = abs(y1).^2 / (n * sum(window.^2));
% f = (0:floor(n/2))/(n*dt);
% P = P(1:floor(n/2)+1);
% P(2:end-1) = 2*P(2:end-1);
% 
% plot(f, P);
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Power Spectrum of Oscillator');
% xlim([0 1]);   

%% Visualize the phase space %%
figure;
plot(x, y);
xlabel('x');
ylabel('y');
title('Phase space of oscillator (5,5)');

%% Visualize the real part %%
% figure;
% plot(t, x);
% xlabel('t');
% ylabel('x');
% title('Re part of oscillator (5,5)');
function Z = initialize_system(m, n, N)
    X = zeros(N, N);
    Y = zeros(N, N);
    for i=1:N
        for j=1:N
            X(i, j) = cos((m*pi*i)/N-1)*(cos((n*pi*j)/N-1));
        end
    end
    Z = [X(:); Y(:)];
end

function [xi_approx, T_approx] = find_upo_candidate(crossing_states, crossing_times)
    min_period = 5.0;
    max_period = 50.0;
    dist_threshold = 2;
    
    best_dist = inf;
    best_pair = [0, 0];
    
    for i = 1:size(crossing_states, 2)
        for j = i+1:size(crossing_states, 2)
            dt = crossing_times(j) - crossing_times(i);
            if dt < min_period || dt > max_period
                continue;
            end
            
            dist = norm(crossing_states(:,i) - crossing_states(:,j));
            if dist < dist_threshold && dist < best_dist
                best_dist = dist;
                best_pair = [i, j];
                T_approx = dt;
            end
        end
    end
    
    if best_dist == inf
        error('No UPO candidates found. Increase simulation time or relax thresholds.');
    end
    
    % Average the two states for better approximation
    xi_approx = (crossing_states(:,best_pair(1)) + crossing_states(:,best_pair(2))) / 2;
end

function [crossing_states, crossing_times] = record_poincare_crossings(N, tspan)
    dt = 0.01;
    a = [ones(1, N^2)/2, zeros(1, N^2)];
    b = 0;
    
    % Initialize chaotic state
    X = 2*rand(N, N) - 1;
    Y = 2*rand(N, N) - 1;
    Z = [X(:); Y(:)];
    
    % Storage
    crossing_states = [];
    crossing_times = [];
    last_crossing = -inf;
    min_interval = 1.0;  % Minimum time between crossings
    
    % Time vector
    t = tspan(1):dt:tspan(2);
    
    for i = 1:length(t)
        % Current section value
        S_prev = dot(a, Z) - b;
        
        % RK4 integration
        k1 = network_dynamics(0, Z, N);
        k2 = network_dynamics(0, Z + dt*k1/2, N);
        k3 = network_dynamics(0, Z + dt*k2/2, N);
        k4 = network_dynamics(0, Z + dt*k3, N);
        Z_next = Z + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        
        S_next = dot(a, Z_next) - b;
        
        % Detect crossing (negative to positive)
        if (S_prev < 0 && S_next >= 0) && (t(i) - last_crossing > min_interval)
            % Linear interpolation to crossing point
            theta = -S_prev / (S_next - S_prev);
            Z_cross = Z + theta*(Z_next - Z);
            
            % Record crossing
            crossing_states(:, end+1) = Z_cross;
            crossing_times(end+1) = t(i) + theta*dt;
            last_crossing = t(i);
        end
        
        Z = Z_next;
    end
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
