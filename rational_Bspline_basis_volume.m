%% Bspline basis volume
clc;clear
num_PP_u = 2;
num_PP_v = 2;
num_PP_w = 2;

mesh_grid_x = 10 * num_PP_u + 1;
mesh_grid_y = 10 * num_PP_v + 1;
mesh_grid_z = 10 * num_PP_w + 1;

xx = linspace(0, 0.99, mesh_grid_x);
yy = linspace(0, 0.99, mesh_grid_y);
zz = linspace(0, 0.99, mesh_grid_z);

CC = zeros(3, mesh_grid_x, mesh_grid_y, mesh_grid_z);

order_u = 1;
order_v = 1;
order_w = 1;

uu = [0, 1];
if (order_u >= 1)
    for ii = 1 : order_u
        uu = [0, uu];
        uu = [uu, 1];
    end
    for ii = 1 : num_PP_u-order_u-1
        uu = [uu(1:order_u+ii), ii / (num_PP_u-order_u) , uu(order_u+ii+1:end)];
    end
end

vv = [0, 1];
if (order_v >= 1)
    for ii = 1 : order_v
        vv = [0, vv];
        vv = [vv, 1];
    end
    for ii = 1 : num_PP_v-order_v-1
        vv = [vv(1:order_v+ii), ii / (num_PP_v-order_v) , vv(order_v+ii+1:end)];
    end
end

ww = [0, 1];
if (order_w >= 1)
    for ii = 1 : order_w
        ww = [0, ww];
        ww = [ww, 1];
    end
    for ii = 1 : num_PP_w-order_w-1
        ww = [ww(1:order_w+ii), ii / (num_PP_w-order_w) , ww(order_w+ii+1:end)];
    end
end

% control points
PP = zeros(3, num_PP_u, num_PP_v, num_PP_w);
PP(:, 1, 1, 1) = [-1; -1; 0];
PP(:, 1, 1, 2) = [-1; 1; 0];
PP(:, 1, 2, 1) = [1; -1; 0];
PP(:, 1, 2, 2) = [1; 1; 0];

PP(:, 2, 1, 1) = [-1; -1; 6];
PP(:, 2, 1, 2) = [-1; 1; 6];
PP(:, 2, 2, 1) = [1; -1; 6];
PP(:, 2, 2, 2) = [1; 1; 6];

% weights
ww_u = zeros(1, num_PP_u);
ww_v = zeros(1, num_PP_v);
ww_w = zeros(1, num_PP_w);
for ii = 1 : num_PP_u
    ww_u(ii) = rand();
end
for ii = 1 : num_PP_v
    ww_v(ii) = rand();
end
for ii = 1 : num_PP_w
    ww_w(ii) = rand();
end
% Bspline
for ll = 1 : mesh_grid_x
    for mm = 1 : mesh_grid_y
        for nn = 1 : mesh_grid_z
            for ii = 0 : num_PP_u-1
                for jj = 0 : num_PP_v-1
                    for kk = 0 : num_PP_w-1
                        CC(:, ll, mm, nn) = CC(:, ll, mm, nn) + Rational_Bspline(ii, order_u, uu, ww_u, xx(ll)) ...
                            .* PP(:, ii+1, jj+1, kk+1) .* Rational_Bspline(jj, order_v, vv, ww_v, yy(mm))...
                            .* Rational_Bspline(kk, order_w, ww, ww_w, zz(nn));
                    end
                end
            end
        end
    end
end

XX = reshape(CC(1,:,:,:),[mesh_grid_x, mesh_grid_y, mesh_grid_z]);
YY = reshape(CC(2,:,:,:),[mesh_grid_x, mesh_grid_y, mesh_grid_z]);
ZZ = reshape(CC(3,:,:,:),[mesh_grid_x, mesh_grid_y, mesh_grid_z]);

for kk = 1 : mesh_grid_z
    mesh(XX(:,:,kk), YY(:,:,kk), ZZ(:,:,kk),LineWidth=1.0);
    hold on
end
hold on
for ii = 1 : num_PP_u
    for jj = 1 : num_PP_v
        for kk = 1 : num_PP_w
            scatter3(PP(1,ii,jj,kk), PP(2,ii,jj,kk), PP(3,ii,jj,kk), 100, 'filled', 'MarkerFaceColor','cyan');
            hold on
        end
    end
end

function RB = Rational_Bspline(ii, pp, uu, ww, xx)
BB = 0.0;
for jj = 0 : length(ww)-1
    BB = BB + Bspline(jj, pp, uu, xx) .* ww(jj+1);
end
if (BB == 0)
    RB = 0.0;
else
    RB = Bspline(ii, pp, uu, xx) .* ww(ii+1) ./ BB;
end
end

function BB = Bspline(ii, pp, uu, xx)
if ( pp == 0 )
    if ( (xx < uu(ii+2)) && (xx >= uu(ii+1)) )
        BB = 1.0;
    else
        BB = 0.0;
    end
else
    a1 = xx - uu(ii+1);
    b1 = uu(ii+pp+1) - uu(ii+1);
    a2 = uu(ii+pp+2) - xx;
    b2 = uu(ii+pp+2) - uu(ii+2);

    if (b1 == 0)
        coe1 = 0;
    else
        coe1 = a1/b1;
    end
    if (b2 == 0)
        coe2 = 0;
    else
        coe2 = a2 / b2;
    end
    BB = coe1 .* Bspline(ii, pp-1, uu, xx) + coe2 .* Bspline(ii+1, pp-1, uu, xx);
end
end
