N = 9;
x = reshape(sol.x(:, 1:N^2), [], N, N);
y = reshape(sol.x(:, N^2 + 1:end), [], N, N);
[X, Y] = meshgrid(1:N, 1:N);
figure;
h1 = surf(X, Y, zeros(N, N), 'EdgeColor', 'none', 'FaceColor', 'blue');
hold on;
h2 = surf(X, Y, zeros(N, N), 'EdgeColor', 'none', 'FaceColor', 'red');
hold off;
zlabel('Re and Img parts of the orbit');
title('Time-Evolving 3D Surface Plot');

shading interp;
view(45, 80); % Set 3D perspective
zlim([-1.5, 1.5]); % Fix z-axis limits for consistency

for i = 1:20:length(sol.t)
    Z_re = squeeze(x(i,:,:));
    Z_img = squeeze(y(i,:,:));
    set(h1, 'ZData', Z_re);
    set(h2, 'ZData', Z_img);
    title(sprintf('Time-Evolving 3D Surface Plot (t = %.1f)', sol.t(i)));
    drawnow;
end


