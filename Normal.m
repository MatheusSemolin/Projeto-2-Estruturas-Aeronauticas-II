function [sec] = Normal(sec)
% Autor: Matheus Semolin Alves Farias
% Normal: Calcula a tensão normal de flexão em cada nó de seções abertas e
% plota o diagrama correspondente.
%
%   A função recebe a struct 'sec' com a geometria, inércias e momentos
%   fletores da seção transversal e retorna a estrutura atualizada.
%
%   Campos obrigatórios na struct de entrada:
%   sec.nodes      - Matriz de nós no formato [ID, x, y, Area_do_Boom].
%                    OBS: A área do boom não é utilizada neste cálculo, pois
%                    os reforçadores já devem estar inclusos nas inércias.
%   sec.panels     - Matriz de painéis no formato [ID_Origem, ID_Destino, Espessura]
%   sec.cargas.Mx         - Momento fletor aplicado em torno do eixo x
%   sec.cargas.My         - Momento fletor aplicado em torno do eixo y
%   sec.propriedades.Ixx        - Momento de inércia em relação ao eixo x
%   sec.propriedades.Iyy        - Momento de inércia em relação ao eixo y
%   sec.propriedades.Ixy        - Produto de inércia
%   sec.propriedades.centroide  - Coordenadas do centróide [xc, yc]
%
%   Saídas:
%   sec.resultados.Sigma - Vetor com os valores da tensão normal em cada nó.


nodes = sec.malha.nodes;
panels = sec.malha.panels;
Mx = sec.cargas.Mx;
My = sec.cargas.My;

N = size(nodes, 1);
P = size(panels, 1);

Ixx = sec.propriedades.Ixx;
Iyy = sec.propriedades.Iyy;
Ixy = sec.propriedades.Ixy;
xc =  sec.propriedades.centroide(1);
yc =  sec.propriedades.centroide(2);

% Coeficientes
Delta = (Ixx*Iyy - Ixy^2);
Kx = (My * Ixx - Mx * Ixy) / Delta;
Ky = (Mx * Iyy - My * Ixy) / Delta;

% Calcula a tensão normal diretamente em cada Nó (coordenadas x, y)
sigma_no = zeros(N, 1);
% fprintf('\n--- TENSÃO NORMAL NOS NÓS (MPa) ---\n');
for i = 1:N
    xn = nodes(i,2) - xc;
    yn = nodes(i,3) - yc;

    % Equação direta da reta de flexão
    sigma_no(i) = Ky * yn + Kx * xn;
    % fprintf('Nó %d: Sigma = %.4f MPa\n', nodes(i,1), sigma_no(i));
end

% Salva na struct
sec.resultados.Sigma_nodes = sigma_no;

% Gráfico do Diagrama de Tensão Normal
figure('Name', 'Diagrama de Tensão Normal', 'Color', 'w');
hold on; grid on; axis equal;
title('Diagrama da Tensão Normal \sigma (MPa)');
xlabel('x (mm)'); ylabel('y (mm)');

% Desenha a seção
for p = 1:P
    n1 = find(nodes(:,1) == panels(p,1));
    n2 = find(nodes(:,1) == panels(p,2));
    plot([nodes(n1,2), nodes(n2,2)], [nodes(n1,3), nodes(n2,3)], 'k', 'LineWidth', 1.5);
end
scatter(nodes(:,2), nodes(:,3), 40, 'k', 'filled');

for i = 1:N
    text(nodes(i,2) + 1.5, nodes(i,3) + 1.5, sprintf('  %d', nodes(i,1)), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'b');
end

% Define a escala
sigma_max = max(abs(sigma_no));
max_dim = max(max(nodes(:,2))-min(nodes(:,2)), max(nodes(:,3))-min(nodes(:,3)));
if sigma_max == 0
    escala = 0.2 * max_dim;
else
    escala = (0.2 * max_dim) / sigma_max;
end

% Plota o diagrama de tensões para cada painel
for p = 1:P
    id_n1 = panels(p,1);
    id_n2 = panels(p,2);

    idx1 = find(nodes(:,1) == id_n1);
    idx2 = find(nodes(:,1) == id_n2);

    x1 = nodes(idx1,2); y1 = nodes(idx1,3);
    x2 = nodes(idx2,2); y2 = nodes(idx2,3);
    dx = x2 - x1; dy = y2 - y1;
    L = sqrt(dx^2 + dy^2);

    % Vetor normal bruto
    nx = -dy / L;
    ny =  dx / L;

    % --- CORREÇÃO: Garante que o vetor normal aponta SEMPRE para fora ---
    x_mid = (x1 + x2) / 2;
    y_mid = (y1 + y2) / 2;
    if (nx * (x_mid - xc) + ny * (y_mid - yc)) < 0
        nx = -nx;
        ny = -ny;
    end

    sig1 = sigma_no(idx1);
    sig2 = sigma_no(idx2);

    x_sig1 = x1 + sig1 * escala * nx;
    y_sig1 = y1 + sig1 * escala * ny;

    x_sig2 = x2 + sig2 * escala * nx;
    y_sig2 = y2 + sig2 * escala * ny;

    % Muda a cor conforme o sinal
    if sig1 * sig2 < 0
        % Encontra a fração 't' onde a tensão é zero
        t = abs(sig1) / (abs(sig1) + abs(sig2));

        x_zero = x1 + t * (x2 - x1);
        y_zero = y1 + t * (y2 - y1);

        % Nó 1 até o Ponto Zero
        sig_media_1 = (sig1 + 0) / 2;
        if sig_media_1 >= 0, cor_1 = 'g'; else, cor_1 = 'r'; end
        plot([x_sig1, x_zero], [y_sig1, y_zero], cor_1, 'LineWidth', 1.5);

        % Ponto Zero até o Nó 2
        sig_media_2 = (0 + sig2) / 2;
        if sig_media_2 >= 0, cor_2 = 'g'; else, cor_2 = 'r'; end
        plot([x_zero, x_sig2], [y_zero, y_sig2], cor_2, 'LineWidth', 1.5);
    else
        % Sem mudança de sinal
        if (sig1 + sig2) / 2 >= 0
            cor = 'g';
        else
            cor = 'r';
        end
        plot([x_sig1, x_sig2], [y_sig1, y_sig2], cor, 'LineWidth', 1.5);
    end
end