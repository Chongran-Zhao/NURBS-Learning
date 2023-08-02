%% Rational Bezier curve
clc;clear
uu = linspace(0, 1, 51);
CC = zeros(3, length(uu));
kk = 3;
PP = zeros(3, kk+1);
PP(:, 1) = [0; 0; 0];
PP(:, 2) = [1; 0; 1];
PP(:, 3) = [0; 1; 1];
PP(:, kk+1) = [1; 1; 1];
ww = [1, 7, 7, 1];

for jj = 1 : 51
    for ii = 1 : kk+1
        CC(:, jj) = CC(:, jj) + Rational_B_basis(ii, kk+1, ww, uu(jj)) .* PP(:, ii);
    end
end

plot3(CC(1,:), CC(2,:), CC(3,:), LineWidth=2.5, Color='r');
hold on
for ii = 1 : kk+1
    scatter3(PP(1,ii), PP(2,ii), PP(3,ii),70, 'filled', 'MarkerFaceColor','cyan');
    hold on
end
legend('k-th order Rational Bezier curve', 'controal point1', 'controal point2', 'controal point3', 'controal point4');

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