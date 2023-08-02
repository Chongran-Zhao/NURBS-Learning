%% Bezier surface
clc;clear
mesh_grid_x = 11;
mesh_grid_y = 11;
xx = linspace(0, 1, mesh_grid_x);
yy = linspace(0, 1, mesh_grid_y);
CC = zeros(3, mesh_grid_x, mesh_grid_y);
order_u = 1;
order_v = 1;
PP = zeros(3, order_u+1, order_v+1);

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        PP(1, ii, jj) = 2 - ii;
        PP(2, ii, jj) = 2 - jj;
        PP(3, ii, jj) = 0;
    end
end
PP(3, 2, 2) = 3;

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        for kk = 1 : mesh_grid_x
            for ll = 1 : mesh_grid_y
                CC(:, kk, ll) = CC(:, kk, ll) + B_basis(ii, order_u+1, xx(kk)) .* PP(:, ii, jj) .* B_basis(jj, order_v+1, yy(ll));
            end
        end
    end
end

XX = reshape(CC(1,:,:),[mesh_grid_x,mesh_grid_y]);
YY = reshape(CC(2,:,:),[mesh_grid_x,mesh_grid_y]);
ZZ = reshape(CC(3,:,:),[mesh_grid_x,mesh_grid_y]);

mesh( XX, YY, ZZ, LineWidth=3.0);
hold on
for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        scatter3(PP(1,ii,jj), PP(2,ii,jj), PP(3,ii,jj), 100, 'filled', 'MarkerFaceColor','cyan');
        hold on
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