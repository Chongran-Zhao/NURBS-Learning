%% Power basis function
% C_0 = a_0;
% C_1 = a_0 + a_1 * u;
% C_2 = a_0 + a_1 * u + a_2 * u^2;
% ......
% C_k = a_0 + a_1 * u + ... + a_k * u^k;
%
% C_k+1 = C_k + a_k+1 * u^k+1;

clc;clear;
mesh_grid_u = 11;
mesh_grid_v = 11;
uu = linspace(0, 1, mesh_grid_u);
vv = linspace(0, 1, mesh_grid_v);

CC = zeros(3, mesh_grid_u, mesh_grid_v);
order_u = 1;
order_v = 1;
PP = zeros(3, order_u+1, order_v+1);

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        PP(:, ii, jj) = 2.0 .* rand(3, 1) + [-1; -1; -1];
        PP(:, ii, jj) = PP(:, ii, jj) ./ norm( reshape(PP(:, ii, jj), [3,1]) );
    end
end


for ii = 1 : mesh_grid_u
    for jj = 1 : mesh_grid_v
        for kk = 1 : order_u+1
            for ll = 1 : order_v+1
                CC(:, ii, jj) = CC(:, ii, jj) + PP(:, kk, ll) .* (uu(ii).^(kk-1)) .* (vv(jj).^(ll-1));
            end
        end
    end
end

XX = reshape(CC(1,:,:),[mesh_grid_u,mesh_grid_v]);
YY = reshape(CC(2,:,:),[mesh_grid_u,mesh_grid_v]);
ZZ = reshape(CC(3,:,:),[mesh_grid_u,mesh_grid_v]);

mesh( XX, YY, ZZ, LineWidth=3.0);
hold on
for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        scatter3(PP(1,ii,jj), PP(2,ii,jj), PP(3,ii,jj), 70, 'filled', 'MarkerFaceColor','cyan');
        hold on
    end
end