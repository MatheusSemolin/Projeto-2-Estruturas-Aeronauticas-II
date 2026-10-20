function [sec] = cisalhamento_aberto_teste(sec)

%% ============================================================
%  DADOS
% =============================================================

nodes  = sec.malha.nodes;
panels = sec.malha.panels;

Vx = sec.cargas.Vx;
Vy = sec.cargas.Vy;

N = size(nodes,1);
P = size(panels,1);

% Propriedades REAIS da seção
Ixx = sec.propriedades.Ixx;
Iyy = sec.propriedades.Iyy;
Ixy = sec.propriedades.Ixy;

xc = sec.propriedades.centroide(1);
yc = sec.propriedades.centroide(2);

%% ============================================================
%  COEFICIENTES DE CISALHAMENTO
% =============================================================

Delta = Ixx*Iyy - Ixy^2;

Qx = (Vx*Ixx - Vy*Ixy)/Delta;
Qy = (Vy*Iyy - Vx*Ixy)/Delta;

fprintf('\n');
fprintf('============================================================\n');
fprintf('             CISALHAMENTO - SECAO ABERTA\n');
fprintf('============================================================\n');

fprintf('\nPROPRIEDADES REAIS DA SECAO\n');
fprintf('Ixx = %.12e mm^4\n',Ixx);
fprintf('Iyy = %.12e mm^4\n',Iyy);
fprintf('Ixy = %.12e mm^4\n',Ixy);
fprintf('xc  = %.12f mm\n',xc);
fprintf('yc  = %.12f mm\n',yc);

fprintf('\nCARREGAMENTO\n');
fprintf('Vx = %.12f N\n',Vx);
fprintf('Vy = %.12f N\n',Vy);

fprintf('\nCOEFICIENTES\n');
fprintf('Qx = %.12e 1/mm^3\n',Qx);
fprintf('Qy = %.12e 1/mm^3\n',Qy);

%% ============================================================
%  VARIAVEIS
% =============================================================

syms s real

q_sym   = sym(zeros(P,1));
q_final = zeros(P,1);

% Primeiros momentos da representacao idealizada
Sx_boom = 0;
Sy_boom = 0;

Sx_painel = 0;
Sy_painel = 0;

% Area da representacao idealizada
A_boom = 0;
A_painel = 0;

%% ============================================================
%  PRIMEIROS MOMENTOS DOS BOOMS
% =============================================================

for i = 1:N

    x = nodes(i,2) - xc;
    y = nodes(i,3) - yc;

    A = nodes(i,4);

    A_boom = A_boom + A;

    Sx_boom = Sx_boom + A*x;
    Sy_boom = Sy_boom + A*y;

end

%% ============================================================
%  FLUXO BASICO q(s)
% =============================================================

for p = 1:P

    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));

    t = panels(p,3);

    % Coordenadas relativas ao centroide REAL
    x1 = nodes(n1,2) - xc;
    y1 = nodes(n1,3) - yc;

    x2 = nodes(n2,2) - xc;
    y2 = nodes(n2,3) - yc;

    % Comprimento do painel
    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    panels(p,4) = L;

    % Coordenadas ao longo do painel
    x_s = ((x2-x1)/L)*s + x1;
    y_s = ((y2-y1)/L)*s + y1;

    % ========================================================
    % Fluxo basico
    % ========================================================

    q_sym(p) = ...
        -t*Qx*int(x_s,s) ...
        -t*Qy*int(y_s,s);

    % Valor do incremento no final do painel
    q_final(p) = double(subs(q_sym(p),s,L));

    % ========================================================
    % Area do painel
    % ========================================================

    A_painel = t*L;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    A_painel = max(A_painel,0);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = real(A_painel);

    A_painel = double(A_painel);

    A_painel = t*L;

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    A_painel = t*L;

    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    % Soma da area dos paineis
    A_painel_total = t*L;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    % Acumula area
    A_painel = t*L;
    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    A_painel = real(A_painel);

    A_painel = t*L;

    A_painel_total = A_painel;

    A_painel = double(A_painel);

    % Acumulacao correta
    A_painel = t*L;
    A_painel_total = A_painel;

    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = real(A_painel);

    A_painel_total = A_painel;

    % Soma global
    A_painel = t*L;
    A_painel = double(A_painel);

    A_painel_total = A_painel;

    A_painel = A_painel_total;

    A_painel = double(A_painel);

    % ========================================================
    % Primeiro momento do painel
    %
    % integral x ds = L*(x1+x2)/2
    % integral y ds = L*(y1+y2)/2
    % ========================================================

    Sx_painel = Sx_painel + ...
        t*L*(x1+x2)/2;

    Sy_painel = Sy_painel + ...
        t*L*(y1+y2)/2;

end

%% ============================================================
%  AREA TOTAL DOS PAINEIS
% =============================================================

A_painel = 0;

for p = 1:P

    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));

    t = panels(p,3);

    x1 = nodes(n1,2);
    y1 = nodes(n1,3);

    x2 = nodes(n2,2);
    y2 = nodes(n2,3);

    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    A_painel = A_painel + t*L;

end

%% ============================================================
%  AREA E CENTROIDE DA REPRESENTACAO IDEALIZADA
% =============================================================

A_ideal = A_boom + A_painel;

if A_ideal ~= 0

    xc_ideal = ...
        (A_boom*xc + ...
        A_painel*xc + ...
        Sx_boom + Sx_painel)/A_ideal + xc;

    yc_ideal = ...
        (A_boom*yc + ...
        A_painel*yc + ...
        Sy_boom + Sy_painel)/A_ideal + yc;

else

    xc_ideal = NaN;
    yc_ideal = NaN;

end

%% ============================================================
%  PRIMEIROS MOMENTOS TOTAIS EM RELACAO AO CENTROIDE REAL
% =============================================================

Sx_total = Sx_boom + Sx_painel;
Sy_total = Sy_boom + Sy_painel;

%% ============================================================
%  INCREMENTOS DE FLUXO NOS BOOMS
% =============================================================

delta_q_boom = zeros(N,1);

for i = 1:N

    xref = nodes(i,2) - xc;
    yref = nodes(i,3) - yc;

    A = nodes(i,4);

    delta_q_boom(i) = ...
        -A*(Qx*xref + Qy*yref);

end

%% ============================================================
%  MATRIZ DE INCIDENCIA
% =============================================================

C = zeros(N,P);

for p = 1:P

    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));

    C(n1,p) = 1;
    C(n2,p) = -1;

end

%% ============================================================
%  VETOR R
% =============================================================

R = zeros(N,1);

for i = 1:N

    soma_q_final = 0;

    % Paineis que terminam no no
    in_panels = find(C(i,:) == -1);

    for k = 1:length(in_panels)

        p = in_panels(k);

        soma_q_final = ...
            soma_q_final + q_final(p);

    end

    R(i) = ...
        delta_q_boom(i) + soma_q_final;

end

%% ============================================================
%  DIAGNOSTICO DOS PRIMEIROS MOMENTOS
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             GEOMETRIA IDEALIZADA\n');
fprintf('============================================================\n');

fprintf('\nArea dos booms       = %.12e mm^2\n',A_boom);
fprintf('Area dos paineis     = %.12e mm^2\n',A_painel);
fprintf('Area idealizada      = %.12e mm^2\n',A_ideal);

fprintf('\nPrimeiros momentos em relacao ao centroide REAL:\n');

fprintf('Sx_boom              = %+ .12e mm^3\n',double(Sx_boom));
fprintf('Sy_boom              = %+ .12e mm^3\n',double(Sy_boom));

fprintf('Sx_painel            = %+ .12e mm^3\n',double(Sx_painel));
fprintf('Sy_painel            = %+ .12e mm^3\n',double(Sy_painel));

fprintf('\nSx_total             = %+ .12e mm^3\n',double(Sx_total));
fprintf('Sy_total             = %+ .12e mm^3\n',double(Sy_total));

fprintf('\nCentroide REAL:\n');
fprintf('xc_real              = %.12f mm\n',xc);
fprintf('yc_real              = %.12f mm\n',yc);

fprintf('\nCentroide IDEALIZADO:\n');
fprintf('xc_ideal             = %.12f mm\n',xc_ideal);
fprintf('yc_ideal             = %.12f mm\n',yc_ideal);

fprintf('\nDiferenca de centroide:\n');
fprintf('dx                   = %+ .12e mm\n',xc_ideal-xc);
fprintf('dy                   = %+ .12e mm\n',yc_ideal-yc);

%% ============================================================
%  BALANCO GLOBAL DO FLUXO
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             BALANCO GLOBAL DE R\n');
fprintf('============================================================\n');

fprintf('\nSoma delta_q_boom    = %+ .12e\n',sum(delta_q_boom));
fprintf('Soma q_final         = %+ .12e\n',sum(q_final));
fprintf('Soma R               = %+ .12e\n',sum(R));

fprintf('\n');

if abs(sum(R)) < 1e-10
    fprintf('>>> BALANCO GLOBAL FECHADO.\n');
else
    fprintf('>>> BALANCO GLOBAL NAO FECHADO.\n');
    fprintf('>>> Nao existe solucao que satisfaca simultaneamente\n');
    fprintf('>>> as 8 equacoes nodais com este modelo.\n');
end

%% ============================================================
%  ESCOLHA DO NO DE REFERENCIA
% =============================================================

grau = sum(C ~= 0,2);

sem_reforco = (nodes(:,4) == 0);

candidatos = find(grau > 1 & sem_reforco);

if isempty(candidatos)

    error('Nao foi encontrado um no adequado para eliminar.');

end

no_descartar = candidatos(1);

%% ============================================================
%  SOLUCAO PARA q0
% =============================================================

idx_manter = true(N,1);

idx_manter(no_descartar) = false;

C_reduzida = C(idx_manter,:);
R_reduzido = R(idx_manter);

q0 = C_reduzida \ R_reduzido;

%% ============================================================
%  FLUXO TOTAL
% =============================================================

q_total_sym = sym(zeros(P,1));

for p = 1:P

    q_total_sym(p) = ...
        q_sym(p) + q0(p);

end

%% ============================================================
%  ERROS NODAIS
% =============================================================

erro_nodal = C*q0 - R;

fprintf('\n');
fprintf('============================================================\n');
fprintf('             EQUILIBRIO NODAL\n');
fprintf('============================================================\n');

fprintf('\nNo de referencia descartado: %d\n', ...
    nodes(no_descartar,1));

fprintf('\n');

for i = 1:N

    fprintf('No %d:  erro = %+ .12e\n', ...
        nodes(i,1),erro_nodal(i));

end

fprintf('\nNorma do erro nodal = %.12e\n', ...
    norm(erro_nodal));

%% ============================================================
%  q0 DE CADA PAINEL
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             q0 POR PAINEL\n');
fprintf('============================================================\n');

for p = 1:P

    fprintf('Painel %d (%d -> %d): q0 = %+ .12e\n', ...
        p, ...
        panels(p,1), ...
        panels(p,2), ...
        q0(p));

end

%% ============================================================
%  VERIFICACAO DA BORDA LIVRE
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('             BORDAS LIVRES\n');
fprintf('============================================================\n');

for i = 1:N

    if grau(i) == 1

        p = find(C(i,:) ~= 0);

        if ~isempty(p)

            q_borda = double(subs(q_total_sym(p),s,0));

            fprintf('No %d - painel %d - q(0) = %+ .12e\n', ...
                nodes(i,1),p,q_borda);

            q_borda_final = double( ...
                subs(q_total_sym(p),s,panels(p,4)));

            fprintf('No %d - painel %d - q(L) = %+ .12e\n', ...
                nodes(i,1),p,q_borda_final);

        end

    end

end

%% ============================================================
%  TENSOES DE CISALHAMENTO
% =============================================================

tau_max = zeros(P,1);

for p = 1:P

    qfun = q_total_sym(p);

    L = panels(p,4);

    q0p = q0(p);

    q_inicio = double(subs(qfun,s,0));

    q_fim = double(subs(qfun,s,L));

    tau_inicio = q_inicio/panels(p,3);

    tau_fim = q_fim/panels(p,3);

    tau_max(p) = max(abs([tau_inicio tau_fim]));

end

%% ============================================================
%  RESULTADOS
% =============================================================

sec.malha.panels = panels;

sec.resultados.Qx = Qx;
sec.resultados.Qy = Qy;

sec.resultados.q0 = q0;

sec.resultados.q_basico = q_sym;
sec.resultados.q_total = q_total_sym;

sec.resultados.delta_q_boom = delta_q_boom;

sec.resultados.C = C;
sec.resultados.R = R;

sec.resultados.erro_nodal = erro_nodal;

sec.resultados.Sx_boom = Sx_boom;
sec.resultados.Sy_boom = Sy_boom;

sec.resultados.Sx_painel = Sx_painel;
sec.resultados.Sy_painel = Sy_painel;

sec.resultados.Sx_total = Sx_total;
sec.resultados.Sy_total = Sy_total;

sec.resultados.A_ideal = A_ideal;

sec.resultados.xc_ideal = xc_ideal;
sec.resultados.yc_ideal = yc_ideal;

sec.resultados.tau_max = tau_max;

%% ============================================================
%  FIM
% =============================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('                     FIM\n');
fprintf('============================================================\n\n');

end