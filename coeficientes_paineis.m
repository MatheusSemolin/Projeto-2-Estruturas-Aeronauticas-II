function sec = coeficientes_paineis(sec, criterio_juncao)
% Identifica a condição de contorno de cada painel e obtém Kb, Kc e Ks.
if nargin < 2, criterio_juncao = 'apoiado'; end
nodes  = sec.malha.nodes;
panels = sec.malha.panels;
L      = sec.propriedades.L_viga;   % distância entre as faces apoiadas [mm]
np = size(panels, 1);

P = struct('Painel',{}, 'nos',{}, 't',{}, 'b',{}, 'razao_ab',{}, 'razao_bc',{}, 'borda1',{}, ...
    'borda2',{}, 'curva_Kb',{}, 'Kb',{}, 'aviso_Kb',{}, 'Kc',{}, 'curva_Ks',{}, 'Ks',{});

notas = {};     % avisos detalhados exibidos abaixo da tabela

for p = 1:np
    n1 = panels(p,1);  n2 = panels(p,2);  t = panels(p,3);

    % largura do painel (comprimento da chapa na seção)
    xy1 = nodes(nodes(:,1) == n1, 2:3);
    xy2 = nodes(nodes(:,1) == n2, 2:3);
    b = hypot(xy2(1) - xy1(1), xy2(2) - xy1(2));

    c = sec.resultados.L_compressao_paineis(p);

    % Razões
    razao_ab = L / b;

    borda1 = tipo_borda(n1, p, panels, criterio_juncao);
    borda2 = tipo_borda(n2, p, panels, criterio_juncao);
    [cKb, cKs] = escolher_curvas(borda1, borda2);

    P(p).Painel = p;
    P(p).nos = [n1 n2];  
    P(p).t = t;  
    P(p).b = b;
    P(p).razao_ab = razao_ab;  
    P(p).borda1 = borda1;  
    P(p).borda2 = borda2;

    % CORREÇÃO 1: Evita divisão por zero se a chapa for tracionada
    % (também cobre c inválido ou negativo: sem compressão, Kb não se aplica)
    if ~isfinite(c) || c <= 0
        P(p).razao_bc = NaN;
        P(p).curva_Kb = '-';
        P(p).Kb = NaN;
        P(p).aviso_Kb = 'sem compressão';
    else
        P(p).razao_bc = b / c;
        P(p).curva_Kb = num2str(cKb);
        % Kb_value limita a razão no último ponto de cada curva e informa se limitou
        [P(p).Kb, limitado, x_max] = Kb_value(cKb, P(p).razao_bc);
        if limitado
            P(p).aviso_Kb = 'LIMITADO';
            notas{end+1} = sprintf( ...
                ' Painel %d (nós %d-%d): b/c = %.2f passa do fim da curva Kb_%s (%.2f). Foi utilizado o último valor da curva Kb = %.2f', ...
                p, n1, n2, P(p).razao_bc, P(p).curva_Kb, x_max, P(p).Kb); %#ok<AGROW>
        else
            P(p).aviso_Kb = '-';
        end
    end

    P(p).Kc = kc_value(borda1, borda2);

    if isempty(cKs)
        P(p).curva_Ks = '-';          % borda livre: sem curva de Ks
        P(p).Ks = NaN;
    else
        P(p).curva_Ks = cKs;
        P(p).Ks = Ks_value(cKs, razao_ab);
    end
end

T = struct2table(P(:));
T = convertvars(T, {'borda1','borda2','curva_Kb','aviso_Kb','curva_Ks'}, 'string');

fprintf('\n Condições de contorno e coeficientes por painel \n\n');
disp(T);
if ~isempty(notas)
    fprintf(' Avisos do Kb:\n');
    fprintf('%s\n', notas{:});
end

% Salva em sec apenas os coeficientes K (um valor por painel, na ordem de sec.malha.panels)
sec.coef_K.Kb = [P.Kb]';
sec.coef_K.Kc = [P.Kc]';       
sec.coef_K.Ks = [P.Ks]';
end

%% ---------- funções auxiliares ----------

function tipo = tipo_borda(no, p, panels, criterio)
    usos = nnz(panels(:,1) == no) + nnz(panels(:,2) == no);
    if usos == 1
        tipo = 'livre';
        return
    end
    outros = setdiff(find(any(panels(:,1:2) == no, 2)), p);   % painéis vizinhos no nó
    switch criterio
        case 'apoiado'
            tipo = 'apoiada';
        case 'engastado'
            tipo = 'engastada';
        case 'espessura'
            if all(panels(outros,3) > panels(p,3))
                tipo = 'engastada';
            else
                tipo = 'apoiada';
            end
        otherwise
            error('Critério de junção desconhecido: "%s".', criterio);
    end
end

function [cKb, cKs] = escolher_curvas(b1, b2)
% Lados de carga sempre apoiados. Escolhe as curvas pelas bordas longitudinais.
%   Kb: 1 = apoiada/apoiada | 2 = engastada/apoiada
%       3 = livre/engastada | 4 = livre/apoiada
%   Ks: 'B' = apoiadas      | 'C' = engastadas   (sem curva se houver borda livre)
    bordas = {b1, b2};
    nl = sum(strcmp(bordas, 'livre'));
    ne = sum(strcmp(bordas, 'engastada'));

    if nl == 2
        error('Painel com as duas bordas longitudinais livres.');
    elseif nl == 1
        cKs = '';
        if ne == 1, cKb = 3; else, cKb = 4; end
    else
        if ne == 0
            cKb = 1;  cKs = 'B';
        elseif ne == 1
            cKb = 2;  cKs = 'B';     % misto: Ks conservador (apoiado)
        else
            cKb = 2;  cKs = 'C';     % sem curva Kb para duas engastadas: usa a 2
        end
    end
end

function Kc = kc_value(b1, b2)
% Valores da tabela 7 e figura 14 do NACA 3781 
% Lados de carga sempre apoiados. 
%   Kc: E = apoiada_apoiada 
%       F = apoiada_engastada
%       B = engastada_engastada
%       F = apoiada_engastada

    tabela_Kc = struct( ...
        'apoiada_apoiada',     4.00, ...
        'apoiada_engastada',   1.28, ...
        'engastada_engastada', 6.98, ...
        'apoiada_livre',       0.43 );
        
    Kc = tabela_Kc.(strjoin(sort({b1, b2}), '_'));
end