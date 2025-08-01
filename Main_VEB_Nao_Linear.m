%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ANÁLISE NÃO LINEAR DE VIGAS DE EULER-BERNOULLI VIA MÉTODO DOS ELEMENTOS FINITOS
%
%   AUTOR: FÁBIO SANTOS
%   DATA: 30/09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parâmetros de limpeza da tela, do Workspace e de figuras, configuração do formato de saída
close all
clear
clc

%% =============================================================================================================
%   PRIMEIRA PARTE: LEITURA DOS NÓS, DA(S) VIGA(S), MATERIAL(AIS), SEÇÃO(ÕES), CARREGAMENTO(S)
%-------------------------------------------------------------------------------------------------------------
LINEAR = true;
% Arquivo
nome = 'Vigas.xlsx';
nosEst = readmatrix(nome, 'Sheet', 'Nós');
vigas = table2cell(readtable(nome, 'Sheet', 'Elementos (Vigas)','VariableNamingRule','preserve'));
materiais = table2cell(readtable(nome, 'Sheet', 'Materiais','VariableNamingRule','preserve'));
secoes = table2cell(readtable(nome, 'Sheet', 'Seções','VariableNamingRule','preserve'));
cargas = readmatrix(nome, 'Sheet', 'Forças');

%% =============================================================================================================
%   SEGUNDA PARTE: DEFINIÇÃO DOS PARÂMETROS DO MÉTODO DE NEWTON-RAPHSON, DAS CARGAS E DA QUADRATURA DE GAUSS
%-------------------------------------------------------------------------------------------------------------
newton = readmatrix(nome, 'Sheet', 'Parâmetros Newton');

es = newton(1);     % Tolerância do método de Newton-Raphson
maxIt = newton(2);  % Número máximo de iterações do método de Newton-Raphson
NLS = newton(3);    % Número de divisões do carregamento (NLS - number load step)
MND = newton(4);    % Número máximo de decrementos de carga (MND - Max number of Decrements)
LDR = newton(5);    % Taxa de decrescimento a cada decremento de carga (em %) (LDR - Load Decrement Rate)
NLGP = newton(6);   % Número de pontos de Gauss para a parcela Linear (NLGP - Number of Linear Gauss Points)
NNLGP = newton(7);  % Número de pontos de Gauss para a parcela não linear (NNLGP - Number of Non Linear Gauss Points)

% Taxa efetiva de decrescimento da carga
ELDR = ((100 - LDR) / 100) .^ (0:MND);       % Effective Load Decrement Rate (number between 0 and 1)

% Armazena os pontos de Gauss a serem calculados
[xLGP, WLGP] = PontosGauss(NLGP);     % xLPG - Linear Gauss Points, WLGP - Weigth Linear Gauss Points
[xNLGP, WNLGP] = PontosGauss(NNLGP);  % xNLGP - Non Linear Gauss Points, WNLGP - Weigth Non Linear Gauss Points

%% =============================================================================================================
%   TERCEIRA PARTE: CRIAÇÃO DO MODELO EM ELEMENTOS FINITOS PARA A ESTRUTURA
%-------------------------------------------------------------------------------------------------------------
% Cria o modelo em elementos finitos
[nosFem, elemFem, iconec] = meshFEM(nosEst, vigas, cargas);

% Determinação do número de graus de liberdade
NDF = size(nosFem, 1) * 3;

% Criação da variável que armazenará os deslocamentos, as deformações e as tensões
% Deslocamentos
u = zeros(NDF, NLS);

% Inicializa as variáveis para armazenar as tensões e deformações
[tensoes, deformacoes] = criaTensaoEDeformacao(elemFem, secoes);

%% =============================================================================================================
%   QUARTA PARTE: EXECUÇÃO DO MÉTODO DE NEWTON PARA CÁLCULO DOS DESLOCAMENTOS DA VIGA
%-------------------------------------------------------------------------------------------------------------
% Definição da carga efetiva que será aplicada sobre a estrutura
delta = (1/NLS:1/NLS:1)';

% Primeiro Laço: Laço que percorrerá todos os passos de carregamento
% Laço que percorrerá os decrementos
controleSolucaoInicialNR = 0;
fprintf('---------------------------------------------------------------------------\n')
fprintf('Começo do processo de análise da Viga Não Linear de Euler-Bernoulli\n')
fprintf('---------------------------------------------------------------------------\n')
for fi = 1:NLS
    fprintf('    Carga %4.d. Razão %4.4f.\n', fi, delta(fi))
    % Ajusta uma variável de controle para o decremento do carregamento
    controleDecremento = false;      % Ativa o decremento
    
    % Variável que armazena o número de decrementos realizados
    kDecremento = 1;
    
    % Matrizes de Rigidez do Elemento para cada elemento para plot da viga
    KauxElem = zeros(6,6,size(elemFem, 1),maxIt);
    fAuxElem = zeros(6,1,size(elemFem, 1),maxIt);

    while true
        % Variável que armazena o número de iterações do Método de Newton-Raphson
        kNewton = 0;

        % Inicialização do erro do Método de Newton-Raphson
        ea = inf;

        % Solução Inicial do Método de Newton-Raphson
        if controleSolucaoInicialNR == 0
            uNR = zeros(NDF, 1);
            uNROld = uNR;
        else
            uNR = u(:,controleSolucaoInicialNR);
            uNROld = uNR;
        end
        
        % Aplicação dos deslocamentos prescritos sobre o vetor solução inicial e cálculo dos graus de
        % liberdade restringidos
        [uNR, glr, glress] = colcacaoCondicoesEssenciais(uNR, nosFem, nosEst);
        fprintf('    ***Começo do Processo Iterativo.***\n')
        while kNewton < maxIt 
            % Controle para testar a qualidade da matriz tangente
            controleMatrizTangente = false;

            fprintf('       Iteração %3.g. ', kNewton)
            % Inicialização da matriz de rigidez global e do vetor de forças global
            K = zeros(NDF, NDF);
            Ktan = K;
            f = zeros(NDF, 1);
            R = f;
            
            % Laço que irá percorrer todos os elementos
            for elem = 1:size(elemFem, 1)
                % Recolhe a conectividade do elemento
                conec = iconec(elem, :)';

                % Recolhe o comprimento do elemento
                Le = elemFem{elem,4};

                % Recolhe qual seção transversal é a do elemento
                secao = elemFem{elem, 11};
                aux = strcmpi(secoes(:,1), secao); 

                % Recolhe a área da seção (Ae) e o momento de Inércia (Ie)
                Ae = secoes{aux, 4};    Ie = secoes{aux, 5};

                % Recolhe qual o material aplicado no elemento
                material = elemFem{elem, 10};
                aux = strcmpi(materiais(:,1), material);

                % Recolhe qual o módulo de elasticidade (Ee)
                Ee = materiais{aux, 2};

                % Definição dos coeficientes das matrizes de rigidez do elemento
                Axx = Ee * Ae;
                Dxx = Ee * Ie;
                
                % Recolhe quais as cargas efetivas do elemento
                fefetiva = [elemFem{elem, 5}, ELDR(kDecremento) * delta(fi) * [elemFem{elem, 6:9}]];

                % Inicialização das matrizes de rigidez do elemento e do vetor de forças equivalentes
                Kelem = zeros(6,6);
                KtanElem = Kelem;
                felem = zeros(6,1);

                %-------------------------------------------------------------------
                % Primeira Integração pela Quadratura de Gauss - PARTE LINEAR
                %-------------------------------------------------------------------
                for igp = 1:NLGP            % [xLGP, WLGP]
                    % Cálculo das funções de forma N e suas derivadas primeira e segunda (dNdx d2Ndx)
                    [N, dNdxi, d2Ndxi] = shapeFunctions(xLGP(igp), Le);

                    % Cálculo da matriz de rigidez K11
                    Kelem([1,2],[1,2]) = Kelem([1,2],[1,2]) + ...
                                         (2/Le) * Axx * transpose(dNdxi([1,2]))*dNdxi([1,2]) * WLGP(igp);

                    % Cálculo da parte linear da matriz de rigidez K22
                    Kelem(3:6,3:6) = Kelem(3:6,3:6) + ...
                                     (8/Le^3) * Dxx * transpose(d2Ndxi(3:6)) * d2Ndxi(3:6) * WLGP(igp);
                    
                    % Cálculo das forças nodais equivalentes àos carregamentos distribuídos
                    % Verifica se o carregamento é linearmente distribuído
                    if fefetiva(1) == 2
                        % Se não for, é uniformemente distribuído
                        felem(1:2,1) = felem(1:2,1) + (Le/2) * transpose(N(1:2)) * fefetiva(2) * WLGP(igp);
                        felem(3:6,1) = felem(3:6,1) + (Le/2) * transpose(N(3:6)) * fefetiva(3) * WLGP(igp);
                    elseif fefetiva(1) == 4
                        % Sendo linearmente distribuído, calcula a força efetiva aplicada na direção x no
                        % ponto da quadratura de Gauss
                        fxx = (fefetiva(2) + fefetiva(4) + (fefetiva(4) - fefetiva(2)) * xLGP(igp)) / 2;
                        % Calcula a integral da força equivalente
                        felem(1:2,1) = felem(1:2,1) + (Le/2) * transpose(N(1:2)) * fxx * WLGP(igp);
                        % Calcula a força efetiva na direção y no ponto da quadratura de Gauss
                        qyy = (fefetiva(3) + fefetiva(5) + (fefetiva(5) - fefetiva(3)) * xLGP(igp)) / 2;
                        % Calcula a integral da força equivalente
                        felem(3:6,1) = felem(3:6,1) + (Le/2) * transpose(N(3:6)) * qyy * WLGP(igp);
                    end                    
                end

                % Armazenamento da parte linear da matriz de rigidez tangente
                KtanElem(1:2,1:2) = Kelem(1:2,1:2);
                KtanElem(3:6,3:6) = Kelem(3:6,3:6);                
                
                %-------------------------------------------------------------------
                % Segunda Integração pela Quadratura de Gauss - PARTE NÃO LINEAR
                %-------------------------------------------------------------------
                if LINEAR == false
                    for igp = 1:NNLGP       % [xNLGP, WNLGP]
                        % Cálculo das funções de forma N e suas derivadas primeira e segunda (dNdx d2Ndx)
                        [N, dNdxi, d2Ndxi] = shapeFunctions(xNLGP(igp), Le);

                        % Cálculo das aproximações de dwdx e dudx
                        dudxi = dNdxi(1:2) * uNR([1,4], 1);
                        dwdxi = dNdxi(3:6) * uNR([2,3,5,6], 1);

                        % Cálculo da matriz não linear de rigidez K12
                        Kelem(1:2, 3:6) = Kelem(1:2, 3:6) + (2/Le^2) * Axx * transpose(dNdxi(1:2)) * dwdxi * ...
                            dNdxi(3:6) * WNLGP(igp);

                        % Cálculo da matriz não linear de rigidez K21
                        Kelem(3:6, 1:2) = Kelem(3:6, 1:2) + (4/Le^2) * Axx * transpose(dNdxi(3:6)) * dwdxi * ...
                            dNdxi(1:2) * WNLGP(igp);

                        % Cálculo da matriz não linear de rigidez K22
                        auxK22 = (4/Le^3) * Axx * transpose(dNdxi(3:6)) * dwdxi^2 * dNdxi(3:6) * WNLGP(igp);
                        Kelem(3:6, 3:6) = Kelem(3:6, 3:6) + auxK22;

                        % Armazenamento das matrizes de rigidez não linear na matriz de rigidez tangente (T12,
                        % T21)
                        KtanElem(1:2, 3:6) = transpose(Kelem(3:6, 1:2));
                        KtanElem(3:6, 1:2) = Kelem(3:6, 1:2);

                        % Cálculo da parcela não linear da matriz de rigidez tangente
                        KtanElem(3:6,3:6) = KtanElem(3:6,3:6) + auxK22 + (8/Le^3) * Axx * transpose(dNdxi(3:6)) * ...
                            ((Le/2)*dudxi + dwdxi^2) * dNdxi(3:6) * WNLGP(igp);
                    end
                end
                % Reorganização da matriz de rigidez do elemento e da matriz tangente de rigidez
                Kelem = Kelem([1,3,4,2,5,6],[1,3,4,2,5,6]);
                KtanElem = KtanElem([1,3,4,2,5,6],[1,3,4,2,5,6]);
                
                KauxElem(:,:,elem,kNewton+1) = Kelem;

                % Reorganização do vetor de forças nodais equivalentes
                felem = felem([1,3,4,2,5,6], 1);
                fAuxElem(:,1,elem,kNewton+1) = felem;
                % Cálculo do Resíduo do Elemento
                Relem = Kelem * uNR(conec) - felem;

                % Sobreposição na matriz de rigidez global e matriz tangente de rigidez global
                K(conec, conec) = K(conec, conec) + Kelem;
                Ktan(conec, conec) = Ktan(conec, conec) + KtanElem;
                f(conec) = f(conec) + felem;
                R(conec) = R(conec) + Relem;
            end

            %------------------------------------------------------------------------------
            %   APLICAÇÃO DAS CARGAS CONCENTRADAS
            %------------------------------------------------------------------------------
            % Controle para recalcular o resíduo
            controleResiduoConcentrado = false;
            % Primeiro carregamento concentrado
            aux = cargas(:,2) == 1;
            if any(aux ~= 0)
                nosPs = cargas(aux, 3);
                for i = 1:size(nosPs)
                    auxNos = sum(nosFem(:, [2 3]) == nosEst(nosPs(i, :), [2 3]), 2) == 2;
                    nosPs(i) = nosFem(auxNos, 1);
                end
                Ps = cargas(aux, 4);
                Psxx = Ps .* cosd(cargas(aux, 6));
                Psyy = Ps .* sind(cargas(aux, 6));
                idxy = [nosPs * 3 - 2, nosPs * 3 - 1];
                f(idxy(:,1)) = f(idxy(:,1)) + Psxx;
                f(idxy(:,2)) = f(idxy(:,2)) + Psyy;
                controleResiduoConcentrado = true;
            end

            % Cargas de momento concentrado
            aux = cargas(:,2) == 3;
            if any(aux ~= 0)
                nosPs = cargas(aux, 3);
                for i = 1:size(nosPs)
                    auxNos = sum(nosFem(:, [2 3]) == nosEst(nosPs(i, :), [2 3]), 2) == 2;
                    nosPs(i) = nosFem(auxNos, 1);
                end
                Ps = cargas(aux, 4);
                idxy = nosPs * 3;
                f(idxy(:,1)) = f(idxy(:,1)) + Ps;
                controleResiduoConcentrado = true;
            end

            % Verifica se precisa recalcular o resíduo
            if controleResiduoConcentrado
                R = K * uNR - f;
            end

            %------------------------------------------------------------------------------
            %   APLICAÇÃO DAS CONDIÇÕES DE CONTORNO ESSENCIAIS PRESCRITAS
            %------------------------------------------------------------------------------
            % Cálculo do vetor de forças com as forças equivalentes dos deslocamentos prescritos
            if ~isempty(glress)
                f = f - sum(K(:, glress) .* uNR(glress)', 2);
            end

            % Junção dos graus de liberdade restringidos
            glr = sort([glr;glress]);

            % Cálculo dos graus de liberdade sem restrição
            glsr = setxor((1:NDF)', glr);

            %------------------------------------------------------------------------------
            %   CÁLCULO DO RESÍDUO E DO VETOR DE INCREMENTO DE SOLUÇÃO
            %------------------------------------------------------------------------------
            % Cálculo do Resíduo
            %R = K(glsr, glsr)*uNR(glsr) - f(glsr);
            
            % Cálculo de uma métrica para o erro pelo viés do Resíduo
            eResiduo = sqrt(sum(R(glsr).^2) / sum(f(glsr).^2));
            fprintf('Resíduo %2.4e\n', eResiduo)
            if eResiduo < es
                u(:, controleSolucaoInicialNR + 1) = uNR;
                controleSolucaoInicialNR = controleSolucaoInicialNR + 1;
                kDecremento = 1;
                break
            end
            
            if rcond(Ktan(glsr, glsr)) < 1e-14
                controleMatrizTangente = true;
                break
            end

            % Cálculo do vetor de incremento da solução
            deltaU = - Ktan(glsr, glsr) \ R(glsr);

            uNR(glsr) = uNR(glsr) + deltaU;

            % Cálculo do erro pelo viés do comprimento do vetor de incrementos em relação à solução
            eDeltaU = norm(uNR(glsr) - uNROld(glsr)) / norm(uNR(glsr));
            fprintf('                     Delta U %2.4e\n', eDeltaU)
            if eDeltaU < es
                u(:, fi) = uNR;
                controleSolucaoInicialNR = controleSolucaoInicialNR + 1;
                kDecremento = 1;
                break
            end
            uNROld = uNR;
            % Atualiza o número de iterações
            kNewton = kNewton + 1;
        end
        fprintf('    ***Fim do Método Iterativo***\n')
        if kNewton == maxIt || controleMatrizTangente
            if controleMatrizTangente
                fprintf('    A matriz tangente é mal condicionada.\n')
            else
                fprintf('    O número máximo de iterações foi atingido sem convergência.\n')
            end
            fprintf('    Haverá Decremento.\n')
            fprintf('    ====Decremento %3.d====\n', kDecremento)
            kDecremento = kDecremento + 1;
        end
        if kDecremento == (MND + 1)
            error('Programa não foi capaz de convergir para uma solução para certa carga!')
        end
        if kDecremento == 1
            fprintf('    Não Haverá Decremento.\n')
            fprintf('---------------------------------------------------------------------------\n\n')
            break;
        end
    end

    % Plot dos deslocamentos da viga
    if fi == 1
        figure('WindowState', 'maximized', 'Units','normalized')
        hold on
        grid on
        box on
        ax = gca;
        Cor = corForcas(fi);
        ytks = DesenhoViga(ax, u(:, controleSolucaoInicialNR), nosFem, elemFem, iconec, Cor, LINEAR);
    else
        Cor = corForcas(fi);
        ytks = DesenhoViga(ax, u(:, controleSolucaoInicialNR), nosFem, elemFem, iconec, Cor, LINEAR, ytks);
    end
end
%% ------------------------------------------------------------------
%                   Começo do Pós-Processamento
%  ------------------------------------------------------------------
% Plot da deflexão da viga para cada incremento de carga













