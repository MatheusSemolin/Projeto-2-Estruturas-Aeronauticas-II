function [sec] = cisalhamento_aberto(sec)
% Autor: Matheus Semolin Alves Farias
% cisalhamento_aberto: Calcula a tensão de cisalhamento de flexão em seções abertas.
%
%   A função recebe a struct 'sec' com a geometria, inércias e esforços
%   cortantes da seção transversal e retorna a estrutura atualizada.
%
%   Campos obrigatórios na struct de entrada:
%   sec.nodes      - Matriz de nós no formato [ID, x, y, Area_do_Boom]
%   sec.panels     - Matriz de painéis no formato [ID_Origem, ID_Destino, Espessura]
%   sec.cargas.Vx         - Força cortante aplicada na direção x
%   sec.cargas.Vy         - Força cortante aplicada na direção y
%   sec.propriedades.Ixx        - Momento de inércia em relação ao eixo x
%   sec.propriedades.Iyy        - Momento de inércia em relação ao eixo y
%   sec.propriedades.Ixy        - Produto de inércia
%   sec.propriedades.centroide  - Coordenadas do centróide [xc, yc]
%
%   Saídas:
%   sec.resultados.Tau.Fletor     - funções simbólicas de Tau nos painéis
%   sec.resultados.Tau.Fletor_max - Vetor com os valores máximos da tensão de
%                                   cisalhamento para cada painel.

nodes = sec.malha.nodes;
panels = sec.malha.panels;
Vx = sec.cargas.Vx;
Vy = sec.cargas.Vy;

N = size(nodes, 1);
P = size(panels, 1);

Ixx = sec.propriedades.Ixx;
Iyy = sec.propriedades.Iyy;
Ixy = sec.propriedades.Ixy;
xc =  sec.propriedades.centroide(1);
yc =  sec.propriedades.centroide(2);

% Q's
Delta = (Ixx*Iyy - Ixy^2);
Qx = (Vx*Ixx - Vy*Ixy) / Delta;
Qy = (Vy*Iyy - Vx*Ixy) / Delta;

% 3. Integração Simbólica ao Longo das Chapas
syms s real
q_sym = sym(zeros(P, 1));

q_final = zeros(P, 1);

for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    t = panels(p,3);

    %Coloca os nós x1,x2,y1,y2 em relação ao centróide
    x1 = nodes(n1,2) - xc;
    y1 = nodes(n1,3) - yc;
    x2 = nodes(n2,2) - xc;
    y2 = nodes(n2,3) - yc;

    L = sqrt((x2-x1)^2 + (y2-y1)^2);

    panels(p,end+1) = L;

    x_s = ((x2 - x1)/L) * s + x1;
    y_s = ((y2 - y1)/L) * s + y1;

    q_sym(p) = -t * Qx * int(x_s, s) - t * Qy * int(y_s, s);
    q_final(p) = double(subs(q_sym(p), s, L));
end

% 4. (Reforçadores)
delta_q_boom = zeros(N, 1);
for i = 1:N
    xref = nodes(i,2) - xc;
    yref = nodes(i,3) - yc;
    A = nodes(i,4);

    delta_q_boom(i) = - A * (Qx * xref + Qy * yref);
end

% 5. Equilíbrio Nodal
C = zeros(N, P); R = zeros(N, 1);
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    C(n1, p) =  1;
    C(n2, p) =  -1;
end

for i = 1:N
    soma_q_final = 0;
    in_panels = find(C(i,:) == -1);
    for k = 1:length(in_panels)
        soma_q_final = soma_q_final + q_final(in_panels(k));
    end

    R(i) = delta_q_boom(i) + soma_q_final;
end

grau = sum(C ~= 0, 2);              % quantos paineis tocam cada no
sem_reforco = (nodes(:,4) == 0);    % nos sem area de boom
candidatos = find(grau > 1 & sem_reforco);

if isempty(candidatos)
    error('Nenhum no ideal para descartar');
end
no_descartar = candidatos(1);

idx_manter = true(N, 1);
idx_manter(no_descartar) = false;
q0 = C(idx_manter, :) \ R(idx_manter);

% 6. Saída Analítica Finais
q_total_sym = sym(zeros(P, 1));
for p = 1:P
    q_total_sym(p) = q_sym(p) + q0(p);
    % fprintf('q0_%d = %.4f N/mm\n', panels(p,1), q0(p));
end

% fprintf('\n--- EQUAÇÕES SIMBÓLICAS DO FLUXO q(s) ---\n');
for p = 1:P
    eq_formatada = char(vpa(expand(q_total_sym(p)/panels(p,3)), 3));
    % fprintf('Painel %d (Nó %d -> Nó %d): Tau(s) = %s\n', p,panels(p,1), panels(p,2), (eq_formatada));
end

erro_numerico = sum(R);
fprintf("%d",erro_numerico)
%% 8. Diagrama de Tensão Projetada na Direção Normal (Usando quiver)
figure('Name', 'Diagrama de Tensão Cisalhante', 'Color', 'w');
hold on; grid on; axis equal;
title('Diagrama da Tensão de Cisalhamento \tau (MPa)');
xlabel('x (mm)'); ylabel('y (mm)');

% Desenha a seção
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    du = nodes(n2,2) - nodes(n1,2);
    dv = nodes(n2,3) - nodes(n1,3);
    plot([nodes(n1,2), nodes(n2,2)], [nodes(n1,3), nodes(n2,3)], 'k','LineWidth', 1.5);
    quiver(nodes(n1,2), nodes(n1,3), du, dv, 0, 'Color', 'k', 'LineWidth', 2, ...
        'ShowArrowHead', 'on', 'MaxHeadSize', 0.5);
end
% --- PLOTAGEM DOS NÓS E REFORÇADORES ---

idx_sem_ref = nodes(:,4) == 0;
idx_com_ref = nodes(:,4) > 0;

if any(idx_sem_ref)
    scatter(nodes(idx_sem_ref, 2), nodes(idx_sem_ref, 3), 40, 'k', 'filled', 'o');
end
for i = 1:N
    text(nodes(i,2) + 1.5, nodes(i,3) + 1.5, sprintf('  %d', nodes(i,1)), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
end
if any(idx_com_ref)
    % Identifica reforçadores diferentes
    areas_unicas = unique(nodes(idx_com_ref, 4));
    num_areas = length(areas_unicas);

    % Paleta de cores dos reforçadores
    if num_areas == 1
        cores = [0, 0.7, 0];
    else
        cores = rand(num_areas, 3);
    end

    for i = 1:num_areas
        idx_area = nodes(:,4) == areas_unicas(i);

        scatter(nodes(idx_area, 2), nodes(idx_area, 3), 60, cores(i,:), ...
            'filled', '^', 'MarkerEdgeColor', 'k');
    end


end
% Escala
tau_max = zeros(P,1);
for p = 1:P
    t = panels(p,3);
    L = panels(p,4);

    q_func = matlabFunction(q_total_sym(p), 'Vars', s);
    range = linspace(0, L, 1000);

    [~, idmax] = max(abs(q_func(range) / t));
    tau_max(p) = q_func(range(idmax));
end

max_dim = max(max(nodes(:,2))-min(nodes(:,2)), max(nodes(:,3))-min(nodes(:,3)));
if max(abs(tau_max)) == 0
    escala = (0.2 * max_dim);
else
    escala = (0.2 * max_dim) / max(abs(tau_max));
end

% Plota os vetores normais
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    t = panels(p,3);

    x1 = nodes(n1,2); y1 = nodes(n1,3);
    x2 = nodes(n2,2); y2 = nodes(n2,3);
    dx = x2 - x1; dy = y2 - y1;
    L = sqrt(dx^2 + dy^2);

    nx = -dy / L;
    ny =  dx / L;

    % --- CORREÇÃO: Garante que o vetor normal aponta SEMPRE para fora ---
    x_mid = (x1 + x2) / 2;
    y_mid = (y1 + y2) / 2;
    if (nx * (x_mid - xc) + ny * (y_mid - yc)) < 0
        nx = -nx;
        ny = -ny;
    end

    num_pts = 20;
    s_vals = linspace(0, L, num_pts);
    q_vals = double(subs(q_total_sym(p), s, s_vals));
    if isscalar(q_vals)
        q_vals = q_vals * ones(size(s_vals));
    end
    tau_vals = q_vals / t;

    x_base = x1 + (dx/L) * s_vals;
    y_base = y1 + (dy/L) * s_vals;

    % Componentes U e V dos vetores
    % Componentes U e V dos vetores
    U = tau_vals * escala * nx;
    V = tau_vals * escala * ny;

    idx_pos = tau_vals >= 0;
    idx_neg = tau_vals < 0;

    % Plot positivos
    quiver(x_base(idx_pos), y_base(idx_pos), U(idx_pos), V(idx_pos), 0, ...
        'Color', 'b', 'ShowArrowHead', 'off', 'LineWidth', 1.2);

    % Plot negativos
    quiver(x_base(idx_neg), y_base(idx_neg), U(idx_neg), V(idx_neg), 0, ...
        'Color', 'r', 'ShowArrowHead', 'off', 'LineWidth', 1.2);

    % Envelope
    x_env = x_base + U;
    y_env = y_base + V;
    plot(x_env, y_env, 'k', 'LineWidth', 1.5);
end
sec.resultados.Tau_fletor_max = tau_max;
sec.resultados.Tau_fletor_sym = q_total_sym;
end