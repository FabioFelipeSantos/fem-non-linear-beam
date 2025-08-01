function [tensoes, deformacoes] = criaTensaoEDeformacao(elemFem, secoes)
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
end