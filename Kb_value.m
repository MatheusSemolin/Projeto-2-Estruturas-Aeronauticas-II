function [Kb, limitado, x_max] = Kb_value(cc, Razao_bc)
% Kb = Kb_value(cc, Razao_bc): valor de Kb na curva cc, razão b/c
% [Kb, limitado, x_max] = Kb_value(cc, Razao_bc): informa também se a razão foi
%   limitada (limitado = true) e qual é o último valor de b/c da curva (x_max).
% Coeficiente Kb obtido por interpolação das curvas dos gráficos.
%
% Entradas:
%   cc
%     1 : Todos os lados apoiados
%     2 : Um engastado e outro apoiados. Lados com carga apoiados.
%     3 : Um livre outro engastado. Lados com carga apoiados.
%     4 : Um livre outro apoiado. Lados com carga apoiados.
%
% Para Razao_bc acima do último ponto de cada curva, Kb assume o último valor
% da curva de ajuste (limite conservador). Se a função for chamada com uma só
% saída, emite um aviso; com mais de uma, só informa em "limitado".
% KB_VALUE(): plota todas as curvas ajustadas.
if nargin == 0
    plotar_todas_curvas();
    return
end
nome_arquivo = fullfile(pasta_dados(), sprintf('Kb_%d.txt', cc));
if ~isfile(nome_arquivo)
    ids = curvas_disponiveis();
    error('Curva Kb_%d não encontrada. Curvas disponíveis: %s', cc, mat2str(ids));
end
s = obter_ajuste(nome_arquivo);
x_max = s.x_max;
limitado = any(Razao_bc(:) > x_max);
Kb = s.f(Razao_bc);
if limitado && nargout < 2
    warning('Kb:foraDoIntervalo', ...
        'Kb_%d: razão b/c = %.2f acima do máximo da curva (%.2f). Usando o último valor do ajuste (limite conservador).', ...
        cc, max(Razao_bc(:)), x_max);
end
end
%% ---------- funções auxiliares ----------
function pasta = pasta_dados()
% pasta 'Kb' ao lado deste arquivo, independente de quem chamou
pasta = fullfile(fileparts(mfilename('fullpath')), 'Kb');
end
function ids = curvas_disponiveis()
lista = dir(fullfile(pasta_dados(), 'Kb_*.txt'));
ids = [];
for i = 1:numel(lista)
    tok = regexp(lista(i).name, '^Kb_(\d+)\.txt$', 'tokens', 'once');
    if ~isempty(tok), ids(end+1) = str2double(tok{1}); end %#ok<AGROW>
end
ids = sort(ids);
end
function s = obter_ajuste(arquivo)
% guarda o ajuste na memória para não refazer o fit a cada chamada
persistent cache
if isempty(cache), cache = containers.Map(); end
if ~isKey(cache, arquivo)
    [f, ~, x_ajuste, ~] = gerar_ajuste(readmatrix(arquivo));
    cache(arquivo) = struct('f', f, 'x_max', max(x_ajuste));
end
s = cache(arquivo);
end
function [f, ajuste, x_ajuste, y_ajuste] = gerar_ajuste(dados)
if size(dados,2) ~= 2
    error('A matriz dados deve possuir exatamente duas colunas: [x y].');
end
x = dados(:,1);  y = dados(:,2);
validos = isfinite(x) & isfinite(y);
x = x(validos);  y = y(validos);
if numel(x) < 4
    error('São necessários pelo menos 4 pontos para realizar o ajuste.');
end
[x, ordem] = sort(x);  y = y(ordem);
ajuste = fit(x, y, 'smoothingspline', 'SmoothingParam', 0.99);

x_max = max(x);
x_ajuste = linspace(min(x), x_max, 500);
y_ajuste = feval(ajuste, x_ajuste);

% acima de x_max: usa o último valor da curva de ajuste
f = @(x_consulta) avaliar_com_assintota(ajuste, x_max, y_ajuste(end), x_consulta);
end
function y = avaliar_com_assintota(ajuste, x_max, y_fim, xq)
y = feval(ajuste, min(xq, x_max));
y(xq > x_max) = y_fim;
y(isnan(xq)) = NaN;      % min(NaN, x_max) devolveria x_max; mantém NaN como NaN
end
function plotar_todas_curvas()
ids = curvas_disponiveis();
figure('Color','w'); hold on;
cores = lines(numel(ids));
for i = 1:numel(ids)
    arquivo = fullfile(pasta_dados(), sprintf('Kb_%d.txt', ids(i)));
    dados = readmatrix(arquivo);
    [~, ~, x_ajuste, y_ajuste] = gerar_ajuste(dados);
    plot(dados(:,1), dados(:,2), 'o', 'Color', cores(i,:), ...
        'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(x_ajuste, y_ajuste, '-', 'Color', cores(i,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Ajuste: Kb\\_%d', ids(i)));
end
grid on; grid minor; box on;
xlabel('\alpha = b/c');  ylabel('K_b');
legend('Location','best');
title('Comparação de Curvas Ajustadas');
set(gca, 'FontSize', 11, 'LineWidth', 1.1);
end