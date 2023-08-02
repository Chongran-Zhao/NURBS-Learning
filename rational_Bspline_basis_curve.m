clc;clear;
for pp = 1 : 6
    num_PP = 10;
    uu = [0, 1];
    if (pp >= 1)
        for ii = 1 : pp
            uu = [0, uu];
            uu = [uu, 1];
        end
        for ii = 1 : num_PP-pp-1
            uu = [uu(1:pp+ii), ii * uu(end)/ (num_PP-pp) , uu(pp+ii+1:end)];
        end
    end

    PP = zeros(3, num_PP);
    PP(:, 1) = [0; 0; 1];
    PP(:, 2) = [2; 2; 1];
    PP(:, 3) = [3; 0; 1];
    PP(:, 4) = [5; 5; 1];
    PP(:, 5) = [7; 6; 1];
    PP(:, 6) = [10; -6; 1];
    PP(:, 7) = [13; 4; 1];
    PP(:, 8) = [14; -2; 1];
    PP(:, 9) = [18; 4; 1];
    PP(:, 10) = [20; 2; 1];

    weight = zeros(1, num_PP);
    for ii = 1 : num_PP
        weight(ii) = rand();
    end

    grid_mesh_x = 100;
    xx = linspace(0, 0.999, grid_mesh_x);
    CC = zeros(3, grid_mesh_x);

    for ii = 1 : grid_mesh_x
        for jj = 0 : num_PP-1
            CC(:, ii) = CC(:, ii) + Rational_Bspline(jj, pp, uu, weight, xx(ii)) .* PP(:, jj+1);
        end
    end

    plot3(CC(1,:), CC(2,:), CC(3,:), LineWidth=2.5, SeriesIndex=pp);
    hold on
    plot3(PP(1,:), PP(2,:), PP(3,:), LineWidth=0.5, LineStyle="--", Color='k');
    grid on
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

    % if (a1 ~= 0) && (b1 == 0)
    %     coe1 = 0.0;
    % elseif (a1 == 0) && (b1 == 0)
    %     coe1 = 1.0;
    % else
    %     coe1 = a1 / b1;
    % end
    % 
    % if (a2 ~= 0) && (b2 == 0)
    %     coe2 = 0.0;
    % elseif (a2 == 0) && (b2 == 0)
    %     coe2 = 1.0;
    % else
    %     coe2 = a2 / b2;
    % end

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

