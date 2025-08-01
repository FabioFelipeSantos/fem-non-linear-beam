function [u, gLR, glress] = colcacaoCondicoesEssenciais(u, nosFem, nosEst)
%------------------------------------------------------------------------------
%   COLOCAÇÃO DAS CONDIÇÕES DE CONTORNO NO VETOR DE SOLUÇÕES INICIAIS
%------------------------------------------------------------------------------
% Apoios aplicados na estrutura
% Restrições na direção x aplicadas aos nós da estrutura 
coordRestringidas = nosEst(nosEst(:,4) == 1, [2, 3]);
if ~isempty(coordRestringidas)
    % Recolhe os nós da malha com essas coordenadas
    nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
    for i = 1:size(coordRestringidas, 1)
        nosMalhaRestringidos(i) = find(sum(nosFem(:,[2,3]) == coordRestringidas(i, :), 2) == 2);
    end

    gLR = 3*nosMalhaRestringidos - 2;
end

% Restrições na direção y aplicadas aos nós da estrutura 
coordRestringidas = nosEst(nosEst(:,5) == 1, [2, 3]);
if ~isempty(coordRestringidas)
    % Recolhe os nós da malha com essas coordenadas
    nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
    for i = 1:size(coordRestringidas, 1)
        nosMalhaRestringidos(i) = find(sum(nosFem(:,[2,3]) == coordRestringidas(i, :), 2) == 2);
    end

    gLR = [gLR; 3*nosMalhaRestringidos - 1];
end

% Restrições de momento aplicados aos nós da estrutura 
coordRestringidas = nosEst(nosEst(:,6) == 1, [2, 3]);
if ~isempty(coordRestringidas)
    % Recolhe os nós da malha com essas coordenadas
    nosMalhaRestringidos = zeros(size(coordRestringidas, 1), 1);
    for i = 1:size(coordRestringidas, 1)
        nosMalhaRestringidos(i) = find(sum(nosFem(:,[2,3]) == coordRestringidas(i, :), 2) == 2);
    end
    
    gLR = [gLR; 3*nosMalhaRestringidos];
end

gLR = sort(gLR);

% Condições essenciais diferentes de apoios (deslocamentos prescritos)
glress = [];
% Deslocamentos prescritos na direção x
aux = ~isnan(nosEst(:, 7));
if any(aux == 1)
    aux = nosEst(aux, 1);
    for i = 1:size(aux,1)
        noFemCondEss = nosFem(sum(nosFem(:, [2,3]) == nosEst(aux(i), [2,3]), 2) == 2, 1);
        idx = noFemCondEss * 3 - 2;
        glress = [glress; idx]; %#ok<*AGROW> 
        u(idx, 1) = nosEst(aux(i), 7);
    end   
end

% Deslocamentos prescritos na direção y
aux = ~isnan(nosEst(:, 8));
if any(aux == 1)
    aux = nosEst(aux, 1);
    for i = 1:size(aux,1)
        noFemCondEss = nosFem(sum(nosFem(:, [2,3]) == nosEst(aux(i), [2,3]), 2) == 2, 1);
        idx = noFemCondEss * 3 - 1;
        glress = [glress; idx];
        u(idx, 1) = nosEst(aux(i), 8);
    end
end

% Deslocamentos prescritos na rotação em torno de z
aux = ~isnan(nosEst(:, 9));
if any(aux == 1)
    aux = nosEst(aux, 1);
    for i = 1:size(aux,1)
        noFemCondEss = nosFem(sum(nosFem(:, [2,3]) == nosEst(aux(i), [2,3]), 2) == 2, 1);
        idx = noFemCondEss * 3;
        glress = [glress; idx];
        u(idx, 1) = nosEst(aux(i), 9);
    end
end
end