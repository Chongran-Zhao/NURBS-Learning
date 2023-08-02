%% Bezier volume
clc;clear
mesh_grid_u = 11;
mesh_grid_v = 11;
mesh_grid_w = 11;

uu = linspace(0, 1, mesh_grid_u);
vv = linspace(0, 1, mesh_grid_v);
ww = linspace(0, 1, mesh_grid_w);

CC = zeros(3, mesh_grid_u, mesh_grid_v, mesh_grid_w);

order_u = 1;
order_v = 1;
order_w = 1;

% control points
PP = zeros(3, order_u+1, order_v+1, order_w+1);
PP(:, 1, 1, 1) = [-1; -1; 0];
PP(:, 1, 1, 2) = [-1; 1; 0];
PP(:, 1, 2, 1) = [1; -1; 0];
PP(:, 1, 2, 2) = [1; 1; 0];

PP(:, 2, 1, 1) = [-1; -1; 6];
PP(:, 2, 1, 2) = [-1; 1; 6];
PP(:, 2, 2, 1) = [1; -1; 6];
PP(:, 2, 2, 2) = [1; 1; 6];

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        for kk = 1 : order_w+1
            for ll = 1 : mesh_grid_u
                for mm = 1 : mesh_grid_v
                    for nn = 1 : mesh_grid_w
                         CC(:, ll, mm, nn) = CC(:, ll, mm, nn) + B_basis(ii, order_u+1, uu(ll))...
                             .* PP(:, ii, jj, kk) .* B_basis(jj, order_v+1, vv(mm)) .* B_basis(kk, order_w+1, ww(nn));
                    end
                end
            end
        end
    end
end

XX = reshape(CC(1,:,:),[mesh_grid_u, mesh_grid_v, mesh_grid_w]);
YY = reshape(CC(2,:,:),[mesh_grid_u, mesh_grid_v, mesh_grid_w]);
ZZ = reshape(CC(3,:,:),[mesh_grid_u, mesh_grid_v, mesh_grid_w]);

for kk = 1 : mesh_grid_w
    mesh(XX(:,:,kk), YY(:,:,kk), ZZ(:,:,kk),LineWidth=1.0);
    hold on
end
hold on
for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        for kk = 1 : order_w+1
            scatter3(PP(1,ii,jj,kk), PP(2,ii,jj,kk), PP(3,ii,jj,kk), 100, 'filled', 'MarkerFaceColor','cyan');
            hold on
        end
    end
end

function BB = B_basis(ii, nn, uu)
if ((uu == 0) && (ii == 1)) || ((uu == 1) && (ii == nn))
    BB = 1.0;
else
    BB = factorial(nn) / ( factorial(ii) * (factorial(nn-ii) )) * (uu.^ii) * ( (1-uu) .^ (nn-ii) );
end
end

% function B_u = B_basis_u(nn, ii, uu)
% if (ii == 1)
%     B_u = nn .* ( -B_basis(ii, nn-1, uu) );
% elseif (ii == nn)
%     B_u = nn .* ( B_basis(ii-1, nn-1, uu) );
% else
%     B_u = nn .* ( B_basis(ii-1, nn-1, uu) - B_basis(ii, nn-1, uu) );
% end
% end