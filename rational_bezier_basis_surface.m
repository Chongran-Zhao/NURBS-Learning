%% Rational Bezier surface
clc;clear
mesh_grid_u = 11;
mesh_grid_v = 11;
uu = linspace(0, 1, mesh_grid_u);
vv = linspace(0, 1, mesh_grid_v);
CC = zeros(3, mesh_grid_u, mesh_grid_v);
order_u = 2;
order_v = 2;
PP = zeros(3, order_u+1, order_v+1);
w_u = [1, 3, 6, 6, 3, 1];
w_v = [1, 3, 6, 6, 3, 1];

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        PP(:, ii, jj) = 2.0 .* rand(3, 1) + [-1; -1; -1];
        PP(:, ii, jj) = PP(:, ii, jj) / norm( reshape(PP(:, ii, jj), [3,1]) );
    end
end

for ii = 1 : order_u+1
    for jj = 1 : order_v+1
        for kk = 1 : mesh_grid_u
            for ll = 1 : mesh_grid_v
                CC(:, kk, ll) = CC(:, kk, ll) + Rational_B_basis(ii, order_u+1, w_u, uu(kk)) .* PP(:, ii, jj) .* Rational_B_basis(jj, order_v+1, w_v, vv(ll));
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
        scatter3(PP(1,ii,jj), PP(2,ii,jj), PP(3,ii,jj), 100, 'filled', 'MarkerFaceColor','cyan');
        hold on
    end
end

function RB = Rational_B_basis(ii, nn, ww, uu)
BB = 0.0;
for jj = 1 : nn
    BB = BB + B_basis(jj, nn, uu) .* ww(jj);
end
if (BB == 0)
    RB = 0.0;
else
    RB = B_basis(ii, nn, uu) .* ww(ii) ./ BB;
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