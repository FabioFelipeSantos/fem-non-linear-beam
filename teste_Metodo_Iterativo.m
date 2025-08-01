%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ANÁLISE NÃO LINEAR DE VIGAS DE EULER-BERNOULLI VIA MÉTODO DOS ELEMENTOS FINITOS
%
%   AUTOR: FÁBIO SANTOS
%   DATA: 30/09/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parâmetros de limpeza da tela, do Workspace e de figuras, configuração do formato de saída
format long
close all
clear
clc

%=============================================================================================================
%   PRIMEIRA PARTE: DEFINIÇÃO DOS NÓS, DA(S) VIGA(S), MATERIAL(AIS), SEÇÃO(ÕES), CARREGAMENTO(S)
%-------------------------------------------------------------------------------------------------------------
% Definição dos nós da estrutura
% nos = [Num Nó, Coord X, Coor Y]
nosEst = [1 0 0; 2 50 0];

% Definição das vigas
% vigas = {Num Viga, Nó inicial, Nó final, Num de Elem Finitos, Material, Seção} (deve ser célula)
vigas = {1, 1, 2, 4, 'M1', 'S1'};

% Definição dos materiais das vigas
% materiais = {Nome, E (módulo de Elasticidade), nu (Coeficiente de Poisson)}
materiais = {'M1', 30e6, 0.3};

% Definição da seção transversal
% secoes = {Nome Seção, Altura (h), Y do Centro de Gravidade (yc), Área (A), Módulo de Inércia (I)}
b = 1;  h = 1; A = b * h; I = (1/12)*b*h^3; yc = h/2;
secoes = {'S1', h, yc, A, I};

% De finição do carregamento
% Tipo -> 1 = Concentrado, 2 = Uniformemente Distribuído e 3 = Momento Concentrado
%         4 = Linearmente Distribuído
% Se tipo = 1, 2 ou 3 -> carga = [Tipo, Nó (se tipo = 1 ou 3) ou Viga, Magnitude, nan, Direção com a viga]
% Se tipo = 4 -> carga = [Tipo, Viga, Magnitude 1 (Nó Inicial), Magnitude 2 (Nó Final), Direção com a viga]
% Defina todas as cargas na mesma variável " cargas "
cargas = [2, 1, 10, nan, 90];

%=============================================================================================================
%   SEGUNDA PARTE: DEFINIÇÃO DAS CONDIÇÕES DE CONTORNO ESSENCIAIS
%-------------------------------------------------------------------------------------------------------------
% Condições essenciais (apoios)
% Tipos -> 1 - Móvel, 2 - Fixo, 3 - Engaste
% Forma de Definir: apoios = [Tipo, Nó, Direção do Apoio]
apoios = [1 1 0];

% Outras formas de condições essenciais
% Forma de Definir: condEss = [Nó, Eixo Global (1 - x, 2 - y, 3 - rotacao z), Valor]
condEss = [2, 1, 0;2 3 0];

%=============================================================================================================
%   TERCEIRA PARTE: DEFINIÇÃO DOS PARÂMETROS DO MÉTODO DE NEWTON-RAPHSON, DAS CARGAS E DA QUADRATURA DE GAUSS
%-------------------------------------------------------------------------------------------------------------
% Tolerância do método de Newton-Raphson
es = 1e-3;

% Número máximo de iterações do método de Newton-Raphson
maxIt = 25;

% Número de divisões do carregamento (NLS - number load step)
NLS = 10;

% Número máximo de decrementos de carga (MND - Max number of Decrements)
MND = 10;

% Taxa de decrescimento a cada decremento de carga (em %) (LDR - Load Decrement Rate)
LDR = 50;

% Taxa efetiva de decrescimento da carga
ELDR = ((100 - LDR) / 100) .^ (0:MND);       % Effective Load Decrement Rate (number between 0 and 1)

% Número de pontos de Gauss para a parcela Linear (NLGP - Number of Linear Gauss Points)
NLGP = 2;

% Número de pontos de Gauss para a parcela não linear (NNLGP - Number of Non Linear Gauss Points)
NNLGP = 1;

% Armazena os pontos de Gauss a serem calculados
[xLGP, WLGP] = PontosGauss(NLGP);     % xLPG - Linear Gauss Points, WLGP - Weigth Linear Gauss Points
[xNLGP, WNLGP] = PontosGauss(NNLGP);  % xNLGP - Non Linear Gauss Points, WNLGP - Weigth Non Linear Gauss Points

%=============================================================================================================
%   QUARTA PARTE: CRIAÇÃO DO MODELO EM ELEMENTOS FINITOS PARA A ESTRUTURA
%-------------------------------------------------------------------------------------------------------------
% Cria o modelo em elementos finitos
[nosFem, elemFem, iconec] = meshFEM(nosEst, vigas, cargas);

% Determinação do número de graus de liberdade
NDF = size(nosFem, 1) * 3;

% Criação da variável que armazenará os deslocamentos, as deformações e as tensões
% Deslocamentos
u = zeros(NDF, NLS);

% Tensões
% Divisões do comprimento do elemento (valores de x ao longo do elemento em que será calculado a tensão)
DEPT = 6;      % Diviões do Elemento em Pontos para a Tensão
xx = [];
for elem = 1:size(elemFem, 1)
    Lelem = elemFem{elem, 4};
    if elem == 1
        xx = (0:Lelem/(DEPT - 1):Lelem)';
    else
        xxInicial = xx(end) + Lelem/(DEPT-1);
        xxFinal = xx(end) + Lelem;
        xx = [xx; (xxInicial:Lelem/(DEPT - 1):xxFinal)']; %#ok<*AGROW> 
    end
end

% Divões ao longo da altura da seção (valores de z ao longo da seção em que serão calculados a tensão)
DZEPT = 5;      % Divisões em Z no Elemento para os Pontos da Tensão
for elem = 1:size(elemFem, 1)
    secao = elemFem{elem, 11};
    aux = strcmpi(secoes(:,1), secao);
    haux = secoes{aux, 2};
    ycaux = secoes{aux, 3};
    zz = (-ycaux:haux/(DZEPT-1):haux-ycaux)';
end

% Cria as variáveis que armazenarão as tensões e as deformações
% tensoes [xx, zz + 1 , elementos]. Na primeira coluna os valores de x em que serão calculados as tensões.
% Nas colunas seguintes da 2 até a zz + 1 os valores da tensão na coordenada z. Cada tensão x por z estará em
% uma página, no total de número de elementos páginas
tensoes = zeros(size(xx, 1), size(zz, 1) + 1, size(elemFem, 1));
tensoes(:,1,:) = repmat(xx, [1,1,size(elemFem, 1)]);

% Deformações. Segue o mesmo processo das tensões
deformacoes = tensoes;

clear Lelem xx xxInicial xxFinal secao aux haux ycaux;      % Limpa algumas variáveis
%=============================================================================================================
%   QUINTA PARTE: EXECUÇÃO DO MÉTODO DE NEWTON PARA CÁLCULO DOS DESLOCAMENTOS DA VIGA
%-------------------------------------------------------------------------------------------------------------

% Definição da carga efetiva que será aplicada sobre a estrutura
delta = (1/NLS:1/NLS:1)';

% Primeiro Laço: Laço que percorrerá todos os passos de carregamento
% Laço que percorrerá os decrementos
controleSolucaoInicialNR = 0;

for fi = 1:NLS
    fprintf('')
    % Ajusta uma variável de controle para o decremento do carregamento
    controleDecremento = false;      % Ativa o decremento
    
    % Variável que armazena o número de decrementos realizados
    kDecremento = 1;

    while true
        % Variável que armazena o número de iterações do Método de Newton-Raphson
        kNewton = 0;

        % Inicialização do erro do Método de Newton-Raphson
        ea = inf;

        % Solução Inicial do Método de Newton-Raphson
        if controleSolucaoInicialNR == 0
            uNR = zeros(NDF, 1);
        else
            uNR = u(:,controleSolucaoInicialNR);
        end
        
        % Aplicação dos deslocamentos prescritos sobre o vetor solução inicial e cálculo dos graus de
        % liberdade restringidos
        [uNR, glr, glress] = colcacaoCondicoesEssenciais(uNR, nosFem, nosEst, condEss, apoios);

        while kNewton < maxIt 
            % Inicialização da matriz de rigidez global e do vetor de forças global
            K = zeros(NDF, NDF);
            Ktan = K;
            f = zeros(NDF, 1);
            
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
                % Reorganização da matriz de rigidez do elemento e da matriz tangente de rigidez
                Kelem = Kelem([1,3,4,2,5,6],[1,3,4,2,5,6]);
                KtanElem = KtanElem([1,3,4,2,5,6],[1,3,4,2,5,6]);

                % Reorganização do vetor de forças nodais equivalentes
                felem = felem([1,3,4,2,5,6], 1);

                % Sobreposição na matriz de rigidez global e matriz tangente de rigidez global
                K(conec, conec) = K(conec, conec) + Kelem;
                Ktan(conec, conec) = Ktan(conec, conec) + KtanElem;
                f(conec) = f(conec) + felem;
            end

            %------------------------------------------------------------------------------
            %   APLICAÇÃO DAS CARGAS CONCENTRADAS
            %------------------------------------------------------------------------------
            % Primeiro carregamento concentrado
            aux = cargas(:,1) == 1;
            if any(aux ~= 0)
                nosPs = cargas(aux, 2);
                for i = 1:size(nosPs)
                    auxNos = sum(nosFem(:, [2 3]) == nosEst(nosPs(i, :), [2 3]), 2) == 2;
                    nosPs(i) = nosFem(auxNos, 1);
                end
                Ps = cargas(aux, 3);
                Psxx = Ps .* cosd(cargas(aux, 5));
                Psyy = Ps .* sind(cargas(aux, 5));
                idxy = [nosPs * 3 - 2, nosPs * 3 - 1];
                f(idxy(:,1)) = f(idxy(:,1)) + Psxx;
                f(idxy(:,2)) = f(idxy(:,2)) + Psyy;
            end

            % Cargas de momento concentrado
            aux = cargas(:,1) == 3;
            if any(aux ~= 0)
                nosPs = cargas(aux, 2);
                for i = 1:size(nosPs)
                    auxNos = sum(nosFem(:, [2 3]) == nosEst(nosPs(i, :), [2 3]), 2) == 2;
                    nosPs(i) = nosFem(auxNos, 1);
                end
                Ps = cargas(aux, 3);
                idxy = nosPs * 3;
                f(idxy(:,1)) = f(idxy(:,1)) + Ps;
            end

            %------------------------------------------------------------------------------
            %   APLICAÇÃO DAS CONDIÇÕES DE CONTORNO ESSENCIAIS PRESCRITAS
            %------------------------------------------------------------------------------
            % Cálculo do vetor de forças com as forças equivalentes dos deslocamentos prescritos
            f = f - sum(K(:, glress) .* uNR(glress)', 2);

            % Junção dos graus de liberdade restringidos
            glr = sort([glr;glress]);

            % Cálculo dos graus de liberdade sem restrição
            glsr = setxor((1:NDF)', glr);

            %------------------------------------------------------------------------------
            %   CÁLCULO DO RESÍDUO E DO VETOR DE INCREMENTO DE SOLUÇÃO
            %------------------------------------------------------------------------------
            % Cálculo do Resíduo
            
            R = K(glsr, glsr)*uNR(glsr) - f(glsr); 
            
            % Cálculo de uma métrica para o erro pelo viés do Resíduo
            eResiduo = sqrt(sum(R.^2) / sum(f(glsr).^2));

            if eResiduo < es
                u(:, controleSolucaoInicialNR + 1) = uNR;
                controleSolucaoInicialNR = controleSolucaoInicialNR + 1;
                break
            end

            uNR(glsr) = K(glsr, glsr) \ f(glsr);

            % Atualiza o número de iterações
            kNewton = kNewton + 1;
        end
        if kNewton == maxIt
            kDecremento = kDecremento + 1;
        end
        if kDecremento == (MND + 1)
            error('Programa não foi capaz de convergir para uma solução para certa carga!')
        end
        if kDecremento == 1
            break;
        end
    end
end















