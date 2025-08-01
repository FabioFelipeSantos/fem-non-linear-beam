function yticks = DesenhoViga(ax, u, nosFem, elemFem, iconec, Cor, LINEAR, ytks)
xi = linspace(-1,1,200)';
%escala = 0.0005;
v = zeros(length(xi), size(elemFem, 1));
for elem = 1:size(elemFem, 1)
    Le = elemFem{elem,4};

    % Armazena todos os GL do elemento
    numGlElemento = iconec(elem, :)';

    % Calcula os pontos ao longo da viga que foram interpolados
    vetor = (1/Le) * (nosFem(elemFem{elem,3}, [2,3]) - nosFem(elemFem{elem,2}, [2,3]));
    pontos = nosFem(elemFem{elem,2}, [2,3]) + vetor .* ((Le/2) * (xi + 1));

    % Calcula as funções de forma do elemento
    [N, ~, ~] = shapeFunctions(xi, Le);
    
    uAux = u(numGlElemento);

    % Calcula os deslocamentos
    ux = N(:,[1,2]) * uAux([1;4]);
    v(:, elem) = N(:,3:6) * uAux([2;3;5;6]);
    
    % Pontos deslocados
    %Taux = [cosSenViga(elem, 1), cosSenViga(elem, 2); -cosSenViga(elem, 2) cosSenViga(elem, 1)];

    pontos = pontos + [ux, v(:, elem)];

    % Faz o desenho da viga
    xElem = [nosFem(elemFem{elem, 2}, 2) nosFem(elemFem{elem, 3}, 2)];
    yElem = [nosFem(elemFem{elem, 2}, 3) nosFem(elemFem{elem, 3}, 3)];
    plot(ax, xElem, yElem,'--b', 'LineWidth', 2.5)
    plot(ax, pontos(:, 1), pontos(:, 2), 'LineWidth', 2, 'Color', Cor)
    
    clear numGlElemento
end
if nargin < 8
    ymin = min(min(v));
    ymax = max(max(v));

    ylim = ax.YLim;

    ticks = linspace(ylim(1), ylim(2), 15);
    tickslabel = arrayfun(@(n){sprintf('%.2e', n)}, linspace(ymin, ymax, 15));
    %tickslabel = {linspace(ymin, ymax, 15)};
    ax.YTick = ticks;
    ax.YTickLabel = tickslabel;

    yticks = [ymin, ymax];
else
    ymin = min(min(min(v)), min(ytks));
    ymax = max(max(max(v)), max(ytks));

    ylim = ax.YLim;

    ticks = linspace(ylim(1), ylim(2), 15);
    tickslabel = arrayfun(@(n){sprintf('%.2e', n)}, linspace(ymin, ymax, 15));
    %tickslabel = {linspace(ymin, ymax, 15)};
    ax.YTick = ticks;
    ax.YTickLabel = tickslabel;
    yticks = [ymin, ymax];
end
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 20;
xlabel('Comprimento $x$ da viga, em metros', 'FontSize',24,'Interpreter','latex','Color','m')
ylabel('Deflex\~{a}o $y$ da viga, em metros', 'FontSize',24,'Interpreter','latex','Color','m')
if LINEAR == 1
    texto = "An\'{a}lise Linear da Viga de Euler-Bernoulli";
else
    texto = "An\'{a}lise N\~{a}o-Linear da Viga de Euler-Bernoulli";
end
title(texto, 'FontSize',28,'Interpreter','latex')
ax.Position = [0.088 0.085 0.9 0.85];
end