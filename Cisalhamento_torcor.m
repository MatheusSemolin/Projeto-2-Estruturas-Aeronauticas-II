function [sec] = Cisalhamento_torcor(sec)
% Autor: Matheus Semolin Alves Farias
% Cisalhamento_torcor: Calcula a tensão de cisalhamento torcional em seções abertas.
%
%   A função recebe a struct 'sec' com a geometria e os esforços da seção
%   transversal e retorna a mesma estrutura atualizada com os resultados.
%
%   Campos obrigatórios na struct de entrada:
%   sec.nodes      - Matriz de nós no formato [ID, x, y, Area_do_Boom]
%   sec.panels     - Matriz de painéis no formato [ID_Origem, ID_Destino, Espessura]
%   sec.cargas.Mz         - Momento torsor aplicado
%   sec.propriedades.centroide  - Coordenadas do centróide (xc, yc)
%   sec.propriedades.centro_de_cisalhamento - Coordenadas (x_se, y_se). 
%                    OBS: Os valores devem estar no sistema global (x,y), 
%                    e NÃO em relação ao centróide.
%
%   Saídas:
%   sec.Tau.torcor - Vetor contendo a tensão de cisalhamento em cada painel (q/t).

nodes = sec.malha.nodes;
panels = sec.malha.panels;
P = length(panels);
N = length(nodes);
A_projetada = 0;
x_se = sec.propriedades.centro_de_cisalhamento(1);
y_se = sec.propriedades.centro_de_cisalhamento(2);
Mz = sec.cargas.Mz;

figure('Name', 'Diagrama de Tensão Cisalhante Tocional', 'Color', 'w');
subplot(1,2,1)
title('Área projetada');
hold on; 
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));

    x_triangulo = [nodes(n1,2), nodes(n2,2),x_se];
    y_triangulo = [nodes(n1,3), nodes(n2,3),y_se];

    Matriz_Triangulo = [x_triangulo',y_triangulo',ones(3,1)];

    A_projetada = A_projetada + 0.5 * det(Matriz_Triangulo);

    x_fechado = [x_triangulo, nodes(n1,2)];
    y_fechado = [y_triangulo, nodes(n1,3)];
    
    fill(x_triangulo, y_triangulo, 'm', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(x_fechado, y_fechado, 'm-', 'LineWidth', 1.5);

    plot([nodes(n1,2), nodes(n2,2)], [nodes(n1,3), nodes(n2,3)], 'k','LineWidth', 1.5);
end
scatter(nodes(:,2), nodes(:,3), 40, 'k', 'filled');
for i = 1:N
    text(nodes(i,2) + 1.5, nodes(i,3) + 1.5, sprintf('  %d', nodes(i,1)), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
end
hold off
q_torcor = Mz/(2*A_projetada);

for p = 1:P
    t = panels(p,3);
    sec.resultados.tau_torcao(p) = q_torcor/t;
end

% =========================================================================
% PLOTAGEM DO DIAGRAMA DE TENSÕES DE CISALHAMENTO (TORÇÃO)
% =========================================================================

% Calcula a escala
max_tau = max(abs(sec.resultados.tau_torcao(1:P)));
max_dim = max(max(nodes(:,2))-min(nodes(:,2)), max(nodes(:,3))-min(nodes(:,3)));

if max_tau == 0
    escala = (0.2 * max_dim); 
else 
    escala = (0.2 * max_dim) / max_tau; 
end

xc = sec.propriedades.centroide(1);
yc = sec.propriedades.centroide(2);



% Desenha a estrutura base
subplot(1,2,2)
hold on; axis equal; grid on;
title('Tensão Cisalhante Torcional')
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    plot([nodes(n1,2), nodes(n2,2)], [nodes(n1,3), nodes(n2,3)], 'k', 'LineWidth', 2);
end
scatter(nodes(:,2), nodes(:,3), 40, 'k', 'filled');

% Plota os vetores e o envelope
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1)); 
    n2 = find(nodes(:,1) == panels(p,2)); 

    x1 = nodes(n1,2); y1 = nodes(n1,3);
    x2 = nodes(n2,2); y2 = nodes(n2,3);

    dx = x2 - x1; dy = y2 - y1;
    L = sqrt(dx^2 + dy^2);

    nx = -dy / L; 
    ny =  dx / L;

    % Garante que aponta para fora
    x_mid = (x1 + x2) / 2;
    y_mid = (y1 + y2) / 2;
    if (nx * (x_mid - xc) + ny * (y_mid - yc)) < 0
        nx = -nx;
        ny = -ny;
    end

    num_pts = 20;
    s_vals = linspace(0, L, num_pts);

    tau_p = sec.resultados.tau_torcao(p);
    tau_vals = tau_p * ones(size(s_vals)); 

    x_base = x1 + (dx/L) * s_vals;
    y_base = y1 + (dy/L) * s_vals;

    U = tau_vals * escala * nx;
    V = tau_vals * escala * ny;

    idx_pos = tau_vals >= 0;
    idx_neg = tau_vals < 0;

    if any(idx_pos)
        quiver(x_base(idx_pos), y_base(idx_pos), U(idx_pos), V(idx_pos), 0, ...
            'Color', 'b', 'ShowArrowHead', 'off', 'LineWidth', 1.2);
    end
    if any(idx_neg)
        quiver(x_base(idx_neg), y_base(idx_neg), U(idx_neg), V(idx_neg), 0, ...
            'Color', 'r', 'ShowArrowHead', 'off', 'LineWidth', 1.2);
    end

    x_env = x_base + U;
    y_env = y_base + V;

    if tau_p >= 0
        plot(x_env, y_env, 'b', 'LineWidth', 1.5);
        plot([x_base(1), x_env(1)], [y_base(1), y_env(1)], 'b', 'LineWidth', 1.5);
        plot([x_base(end), x_env(end)], [y_base(end), y_env(end)], 'b', 'LineWidth', 1.5);
    else
        plot(x_env, y_env, 'r', 'LineWidth', 1.5);
        plot([x_base(1), x_env(1)], [y_base(1), y_env(1)], 'r', 'LineWidth', 1.5);
        plot([x_base(end), x_env(end)], [y_base(end), y_env(end)], 'r', 'LineWidth', 1.5);
    end
end
end