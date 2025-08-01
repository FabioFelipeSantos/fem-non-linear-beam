%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   FUNÇÃO MONTA A MALHA 1D DE ELEMENTOS FINITOS DOS ELEMENTOS DE VIGA
%   
%   ENTRADA: nos - Nos da estrutura
%            vigas - vigas que compõe a estrutura
%
%   SAÍDA: nosFem - os nós da malha em elementos finitos 1D;
%          elemFem - elementos da malha de elementos finitos com seu número, nós inicial e final, material e
%                    seção transversal;
%          iconec - matriz de conectividade de cada elemento.
%
%   AUTOR: FÁBIO SANTOS
%   DATA: 30/09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nosFem, elemFem, iconec] = meshFEM(nos, vigas, cargas)
    % Calcula o total de elementos da estrutura no todo
    numTotalElementos = sum([vigas{:,4}]);
    
    % Calcula o total de nós que terá a estrutura
    numTotalNos = numTotalElementos + 1;

    % Inicialização das matrizes de nos e elementos da malha
    nosFem = nan(numTotalNos, 3);          % [Num Nó, CX, CY]
    elemFem = cell(numTotalElementos, 11);    % [Num Elem, Nó Ini, Nó Final, Comprimento, Tipo carga, Carga 1 
    %                                            x, Carga 1 y, Carga 2 x, Carga 2 y, Material, Seção]
    % Inicialização do vetor de conectividades do elemento
    iconec = zeros(numTotalElementos, 6);

    % Variável de controle da posição do elemento a ser adicionado
    auxPosElem = 1;
    % Variável de controle da posição do(s) nó(s) a ser(em) adicionado(s)
    auxPosNo = 1;

    % Laço que percorre cada uma das vigas entradas
    for v = 1:size(vigas, 1)
        % Recolhe quantos elementos há nessa viga
        numElem = vigas{v, 4};
        % Recolhe as coordenadas do nó inicial da viga
        noInicial = nos(vigas{v,2}, [2,3]);
        % Recolhe as coordenadas do nó final da viga
        noFinal = nos(vigas{v,3}, [2,3]);
        % Cria um vetor que liga o nó inicial ao nó final
        vetor = noFinal - noInicial;
        % Calcula o comprimento da viga
        L = norm(vetor);
        % Normaliza o vetor de direção da viga
        vetor = vetor / L;
        % Calcula o comprimento de cada elemento pertencente à viga
        deltaL = L / numElem;

        % Cria os nós adicionais na viga para acomodar todos os elementos
        auxNos = noInicial + deltaL * (0:numElem)' .* vetor;
        
        % Se for a primeira viga, podemos adicionar todos os nós na matriz de nós da malha
        if v == 1
            % Recolhe todos os índices dos nós adicionados
            idx = (auxPosNo:size(auxNos, 1))';
            % Armazena os índices
            nosFem(idx, 1) = idx;
            % Armazena os nós
            nosFem(idx, [2,3]) = auxNos;
            % Atualiza a posição dos nós
            auxPosNo = idx(end) + 1;
        else
            % Caso contrário, poderá haver algum nó que seja repetido, então não iremos adicioná-lo
            for i = 1:size(auxNos, 1)
                % Variável auxiliar que compara os nós adicionados com o a ser adicionado
                aux = sum(nosFem(:,[2,3]) == auxNos(i,:), 2) == 2;
                if all(aux == 0)        % Se não tiver sido adicionado, entra no condicional
                    % Armazena o número do nó
                    nosFem(auxPosNo, 1) = auxPosNo;
                    % Armazena o nó
                    nosFem(auxPosNo, [2,3]) = auxNos(i,:);
                    % Atualiza a posição
                    auxPosNo = auxPosNo + 1;
                end
            end
        end
        
        % Laço para adicionar os elementos que são criados a partir dos nós entre o nó final e inicial da viga
        for i = 1:(size(auxNos, 1)-1)
            % Recolhe o nó inicial do elemento
            auxNoInicial = auxNos(i, :);
            % Recolhe o índice na matriz de nós da malha
            idxInicial = nosFem(sum(nosFem(:,[2,3]) == auxNoInicial, 2) == 2, 1);
            % Recolhe o nó final do elemento
            auxNoFinal = auxNos(i+1, :);
            % Recolhe o índice na matriz de nós da malha
            idxFinal = nosFem(sum(nosFem(:,[2,3]) == auxNoFinal, 2) == 2, 1);
            % Armazena o número do nó
            elemFem(auxPosElem, 1) = {auxPosElem};
            % Armazena os índices dos nós inicial e final do elemento
            elemFem(auxPosElem, [2, 3]) = num2cell([idxInicial, idxFinal]);
            % Calcula o comprimento do elemento
            elemFem(auxPosElem, 4) = {norm(auxNoFinal - auxNoInicial)};
            % Armazena o material do elemento
            elemFem(auxPosElem, 10) = vigas(v, 5);
            % Armazena a seção do elemento
            elemFem(auxPosElem, 11) = vigas(v, 6);
            % Calcula os graus de liberdade do nó inicial
            glInicial = elemFem{auxPosElem, 2} * 3 - (2:-1:0);
            % Calcula os graus de liberdade do nó final
            glFinal = elemFem{auxPosElem, 3} * 3 - (2:-1:0);
            % Armazena os graus de liberdade na matriz de conectividade
            iconec(auxPosElem, :) = [glInicial, glFinal];

            %-----------------------------------------------------------------------------------------
            %   COLOCA AS CARGAS EM CADA ELEMENTO, SE A VIGA POSSUIR CARGA DISTRIBUÍDA
            %-----------------------------------------------------------------------------------------
            % Verifica se tem cargas distribuídas aplicadas sobre a viga
            aux = sum([cargas(:,3) == v, cargas(:,2) == 2, cargas(:,2) == 4], 2) > 1;
            if any(aux ~= 0)
                % Se tiver, define dois controles: um para carga uniformemente distribuída e outra para
                % linearmente distribuiída
                dist = false;
                lindist = false;
                % Verifica se existe apenas uma carga aplicada sobre a viga
                if sum(aux) == 1
                    % Havendo apenas uma carga, verifica se ela é uniformemente distribuída
                    if cargas(aux, 2) == 2
                        % Se for armazena o valor da magnitude da carga e ativa o controle
                        w = cargas(aux, 4);
                        dist = true;
                    elseif cargas(aux, 2) == 4
                        % Se for linearmente distribuída, verifica se estamos no primeiro elemento da viga
                        if i == 1
                            % Armazena a carga no nó 0
                            w1 = cargas(aux, 4);
                            % Armazena a carga no comprimento do elemento
                            w2 = cargas(aux, 4) + (elemFem{auxPosElem, 4}/L) * (cargas(aux, 5) - cargas(aux, 4));
                            % Ativa o controle
                            lindist = true;
                        else 
                            % Calcula o ponto inicial do elemento atual, que será o ponto final do elemento
                            % anterior
                            ponto1 = sum([elemFem{1:auxPosElem-1, 4}]);
                            % Calcula o ponto final do elemento atual
                            ponto2 = ponto1 + elemFem{auxPosElem, 4};
                            % Armazena a carga no comprimento do elemento anterior
                            w1 = cargas(aux, 4) + (ponto1/L) * (cargas(aux, 5) - cargas(aux, 4));
                            % Armazena a carga no comprimento do elemento atual
                            w2 = cargas(aux, 4) + (ponto2/L) * (cargas(aux, 5) - cargas(aux, 4));
                            % Ativa o controle
                            lindist = true;
                        end
                    end
                    if dist
                        % Armazena as cargas uniformemente distribuídas
                        elemFem(auxPosElem, 5) = {2};
                        elemFem(auxPosElem, 6) = {w * cosd(cargas(aux, 6))};
                        elemFem(auxPosElem, 7) = {w * sind(cargas(aux, 6))};
                        elemFem(auxPosElem, 8) = {nan};
                        elemFem(auxPosElem, 9) = {nan};
                    elseif lindist
                        % Armazena as cargas linearmente distribuídas
                        elemFem(auxPosElem, 5) = {4};
                        elemFem(auxPosElem, 6) = {w1 * cosd(cargas(aux, 6))};
                        elemFem(auxPosElem, 7) = {w1 * sind(cargas(aux, 6))};
                        elemFem(auxPosElem, 8) = {w2 * cosd(cargas(aux, 6))};
                        elemFem(auxPosElem, 9) = {w2 * sind(cargas(aux, 6))};
                    end
                else
                    % Há mais de uma carga no mesmo elemento
                    % Armazena essas cargas
                    cargasaux = cargas(aux, :);
                    % Cria uma variável de controle
                    controleDist = zeros(2, sum(aux));
                    % Cria as somas das cargas nas direções x e y nos nós 1 e 2
                    soma1x = 0;
                    soma1y = 0;
                    soma2x = 0;
                    soma2y = 0;
                    % Laço que irá percorrer todas as cargas no mesmo elemento
                    for k = 1:sum(aux)
                        % Verifica se a carga é uniformemente distribuída
                        if cargasaux(k, 1) == 2
                            % Armazena a carga
                            w = cargasaux(k, 4);
                            % Decompõe a carga nas componentes x e y
                            soma1x = soma1x + w * cosd(cargasaux(k, 6));
                            soma1y = soma1y + w * sind(cargasaux(k, 6));
                            % Ativa o controle
                            controleDist(1, k) = 1;
                        else
                            if i == 1
                                % Armazena a carga no nó 0
                                w1 = cargasaux(k, 4);
                                % Armazena a carga no comprimento do elemento
                                w2 = cargasaux(k, 4) + (elemFem{auxPosElem, 4}/L) * (cargasaux(k, 5) - cargasaux(k, 4));
                            else
                                % Calcula o ponto inicial do elemento atual, que será o ponto final do elemento
                                % anterior
                                ponto1 = elemFem{auxPosElem-1, 4};
                                % Calcula o ponto final do elemento atual
                                ponto2 = elemFem{auxPosElem-1, 4} + elemFem{auxPosElem, 4};
                                % Armazena a carga no comprimento do elemento anterior
                                w1 = cargasaux(k, 4) + (ponto1/L) * (cargasaux(k, 5) - cargasaux(k, 4));
                                % Armazena a carga no comprimento do elemento atual
                                w2 = cargasaux(k, 4) + (ponto2/L) * (cargasaux(k, 5) - cargasaux(k, 4));                                
                            end
                            % Decompõe a carga nas componentes x e y
                            soma1x = soma1x + w1 * cosd(cargasaux(k, 6));
                            soma1y = soma1y + w1 * sind(cargasaux(k, 6));
                            soma2x = soma2x + w2 * cosd(cargasaux(k, 6));
                            soma2y = soma2y + w2 * sind(cargasaux(k, 6));
                            % Ativa o controle
                            controleDist(2, k) = 1;                            
                        end
                    end
                    if all(controleDist(2,:) == 0)
                        % Armazena as cargas uniformemente distribuídas
                        elemFem(auxPosElem, 5) = {2};
                        elemFem(auxPosElem, 6) = {soma1x};
                        elemFem(auxPosElem, 7) = {soma1y};
                        elemFem(auxPosElem, 8) = {nan};
                        elemFem(auxPosElem, 9) = {nan};
                    else
                        % Armazena as cargas linearmente distribuídas
                        elemFem(auxPosElem, 5) = {4};
                        elemFem(auxPosElem, 6) = {soma1x};
                        elemFem(auxPosElem, 7) = {soma1y};
                        elemFem(auxPosElem, 8) = {soma2x};
                        elemFem(auxPosElem, 9) = {soma2y};
                    end
                end
            end    
            % Atualiza a posição do elemento
            auxPosElem = auxPosElem + 1;
        end
    end
end
















