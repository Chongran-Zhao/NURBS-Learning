%% Power basis function
% C_0 = a_0;
% C_1 = a_0 + a_1 * u;
% C_2 = a_0 + a_1 * u + a_2 * u^2;
% ......
% C_k = a_0 + a_1 * u + ... + a_k * u^k;
% 
% C_k+1 = C_k + a_k+1 * u^k+1;

clc;clear;
order = 2;
uu = linspace(0, 1, 11);
AA = zeros(3, order+1);
yy = zeros(3, 11);

for ii = 1 : order+1
    AA(:, ii) = 2.0 .* rand(3, 1) + [-1; -1; -1];
    AA(:, ii) = AA(:, ii) / norm(AA(:, ii));
end
for ii = 1 : 11
    CC = AA(:, 1);
    for jj = 2 : order+1
        CC = CC + AA(:, jj) .* uu(ii) .^ (jj-1);
    end
    yy(:, ii) = CC;
end
CX = AA(:, 2);
for ii = 3 : order+1
    CX = CX + (ii-1) * AA(:, ii) .* uu(end) .^ (ii-2);
end
y_x = CX ./ 20.0;

plot3( yy(1,:), yy(2,:), yy(3,:), LineWidth=3.0, Color='r');
hold on
quiver3( yy(1,end), yy(2,end), yy(3,end), yy(1,end)-yy(1,end-1), yy(2,end)-yy(2,end-1), ...
    yy(3,end)-yy(3,end-1), LineWidth=6.5, Color = 'b');
hold on
quiver3( yy(1,end), yy(2,end), yy(3,end), y_x(1), y_x(2), ...
    y_x(3), LineWidth=3.5, Color = 'g');
legend('Power basis curve', '1st order derivative', 'numerical 1st order derivative');
hold on
for ii = 1 : order+1
    scatter3(AA(1,ii), AA(2,ii), AA(3,ii), 70, 'filled', 'MarkerFaceColor','cyan');
    hold on
end
