Re_Z = sol.z(:, 1:N^2);
Re_W = sol.w(:, 1:N^2);
Img_Z = sol.z(:, N^2+1:end);
Img_W = sol.w(:, N^2+1:end);

x = sum(Re_W.^2 + Img_W.^2, 2)/(N^2);
figure;
plot(x);