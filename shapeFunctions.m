function [N, dNdx, d2Ndx] = shapeFunctions(xi, L)
% Funções de Forma
N(:,1) = (1 - xi) / 2;
N(:,2) = (xi + 1) / 2;
N(:,3) = (xi - 1).^2 .* (xi + 2) / 4;
N(:,4) = L * (xi - 1).^2 .* (xi + 1) / 8;
N(:,5) = -(xi + 1).^2 .* (xi - 2) / 4;
N(:,6) = L * (xi - 1) .* (xi + 1).^2 / 8;

% Derivadas das funções
dNdx(:,1) = -1/2 * ones(length(xi), 1);
dNdx(:,2) = 1/2 * ones(length(xi), 1);
dNdx(:,3) = (3*xi.^2)/4 - 3/4;
dNdx(:,4) = (L*(3*xi.^2 - 2*xi - 1))/8;
dNdx(:,5) = 3/4 - (3*xi.^2)/4;
dNdx(:,6) = (L*(3*xi.^2 + 2*xi - 1))/8;


% Derivadas Segundas
d2Ndx(:,1) = zeros(length(xi), 1);
d2Ndx(:,2) = zeros(length(xi), 1);
d2Ndx(:,3) = (3*xi)/2;
d2Ndx(:,4) = (L*(3*xi - 1))/4;
d2Ndx(:,5) = -(3*xi)/2;
d2Ndx(:,6) = (L*(3*xi + 1))/4;
end