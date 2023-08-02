%% Bezier curve
clc;clear
uu = linspace(0, 1, 51);
CC = zeros(3, length(uu));
kk = 5;
PP = zeros(3, kk+1);
PP(:, 1) = [0; 0; 0];
PP(:, kk+1) = [1; 1; 1];
for ii = 2 : kk
    PP(:, ii) = 2.0 .* rand(3, 1) + [-1; -1; -1];
    PP(:, ii) = PP(:, ii) / norm(PP(:, ii));
end
for ii = 1 : kk+1
    for jj = 1 : 51
        CC(:, jj) = CC(:, jj) + B_basis(ii, kk+1, uu(jj)) .* PP(:, ii);
    end
end
scatter3(PP(1,1), PP(2,1), PP(3,1), 70,'filled', 'MarkerFaceColor','blue');
hold on
plot3(CC(1,:), CC(2,:), CC(3,:), LineWidth=2.5, Color='r');
hold on
scatter3(PP(1,end), PP(2,end), PP(3,end),70, 'filled', 'MarkerFaceColor','cyan');
legend('Fixed start node', 'k-th order Bezier curve', 'Fixed end node');

function BB = B_basis(ii, nn, uu)
BB = factorial(nn) / ( factorial(ii) * (factorial(nn-ii) )) * (uu.^ii) * ( (1-uu) .^ (nn-ii) );
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