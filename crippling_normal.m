function [sec] = crippling_normal(sec)
% CRIPPLING_NORMAL Calculo do Crippling Normal pelo metodo Lockheed
%
%   criterio_juncao - 'apoiado', 'engastado' ou 'espessura' (Padrao: 'apoiado')
% Entradas:
%   sec             - Struct principal da secao contendo malha e propriedades

criterio_juncao = 'apoiado';
KSI2MPA = 6.89476;
% 1. Extração de propriedades e malha
nodes  = sec.malha.nodes;
panels = sec.malha.panels;
np     = size(panels, 1);
N = sec.cargas.N;
A = sec.propriedades.A;
A_elem         = zeros(np, 1);
Fcc_elem       = zeros(np, 1);
tipos_lockheed = zeros(np, 1);

% 2. Loop por painel (Trecho do seu amigo adaptado)
for p = 1:np
    n1 = panels(p,1);  
    n2 = panels(p,2);  
    t  = panels(p,3);

    % Geometria: cálculo da largura plana (b) do painel
    xy1 = nodes(nodes(:,1) == n1, 2:3);
    xy2 = nodes(nodes(:,1) == n2, 2:3);
    b   = hypot(xy2(1) - xy1(1), xy2(2) - xy1(2));

    % Determinação do tipo de borda longitudinal
    borda1 = tipo_borda_crippling(n1, p, panels, criterio_juncao);
    borda2 = tipo_borda_crippling(n2, p, panels, criterio_juncao);

    razao_bt = b/t;
    if strcmp(borda1, 'livre') || strcmp(borda2, 'livre')
        cc = 1;
    else
        cc = 2;
    end
    
    A_elem(p)         = b * t;
    tipos_lockheed(p) = cc;
    Fcc_elem(p)       = interpolar_lockheed(razao_bt,cc, false);
end
    
    A_reforcadores = sum(nodes(:,4));
    A_total = A_reforcadores + sum(A_elem);
% 3. Integração na estrutura global (sec)
sec.crippling.Fcc_elem   = Fcc_elem * KSI2MPA;
sec.crippling.Fcc_perfil = (sum(A_elem.*Fcc_elem)/ A_total)*KSI2MPA; % Fcc total [MPa]
sec.MS.crippling_normal = sec.crippling.Fcc_perfil/(N/A) -1;
end

%% ---------- Função Auxiliar de Borda ----------
function tipo = tipo_borda_crippling(no, p, panels, criterio)
usos = nnz(panels(:,1) == no) + nnz(panels(:,2) == no);
if usos == 1
    tipo = 'livre';
    return;
end
outros = setdiff(find(any(panels(:,1:2) == no, 2)), p);
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