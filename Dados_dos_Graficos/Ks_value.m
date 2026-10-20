function Ks = Ks_value(cc, Razao_bc)
% Coeficiente Ks obtido por interpolação das curvas dos gráficos.
%
% Ks_value(cc, Razao_bc): valor de Ks na curva cc, razão b/c
% Ks_value(): plota todas as curvas ajustadas disponíveis.
%
% Entradas:
%   cc
%     'A' : Todos os lados engastados
%     'B' : Todos os lados apoiados
%     'C' : Dois engastados. Lados de carga apoiados.
%     'D' : Dois apoiados. Lados de carga engastados.

    if nargin == 0
        plotar_todas_curvas();
        return
    end

    [ids, arquivos] = curvas_disponiveis();
    id = char(string(cc));                 % aceita 'A', "A" ou número (3)
    idx = find(strcmpi(ids, id), 1);       % ignora maiúscula/minúscula
    if isempty(idx)
        error('Curva Ks_%s não encontrada. Curvas disponíveis: %s', ...
              id, strjoin(ids, ', '));
    end

    f = obter_ajuste(arquivos{idx});
    Ks = f(Razao_bc);
end

%% ---------- funções auxiliares ----------

function pasta = pasta_dados()
    % pasta 'Ks' ao lado deste arquivo, independente de quem chamou
    pasta = fullfile(fileparts(mfilename('fullpath')), 'Ks');
end

function [ids, arquivos] = curvas_disponiveis()
    pasta = pasta_dados();
    if ~isfolder(pasta)
        error('Pasta de dados não encontrada: %s', pasta);
    end
    lista = dir(fullfile(pasta, 'Ks_*.txt'));
    ids = {};  arquivos = {};
    for i = 1:numel(lista)
        tok = regexp(lista(i).name, '^Ks_(.+)\.txt$', 'tokens', 'once');
        if ~isempty(tok)
            ids{end+1} = tok{1};                                         %#ok<AGROW>
            arquivos{end+1} = fullfile(lista(i).folder, lista(i).name);  %#ok<AGROW>
        end
    end
    if isempty(ids)
        error('Nenhum arquivo Ks_*.txt encontrado em: %s', pasta);
    end
    [ids, ord] = sort(ids);
    arquivos = arquivos(ord);
end

function f = obter_ajuste(arquivo)
    % guarda o ajuste na memória para não refazer o fit a cada chamada
    persistent cache
    if isempty(cache), cache = containers.Map(); end
    if ~isKey(cache, arquivo)
        cache(arquivo) = gerar_ajuste(readmatrix(arquivo));
    end
    f = cache(arquivo);
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

% Assíntota: para x acima do último ponto, mantém o valor em x máximo
x_max = max(x);
f = @(x_consulta) feval(ajuste, min(x_consulta, x_max));

x_ajuste = linspace(min(x), x_max, 500);
y_ajuste = f(x_ajuste);
end

function plotar_todas_curvas()
    [ids, arquivos] = curvas_disponiveis();
    fprintf('Pasta de dados: %s\n', pasta_dados());
    fprintf('Curvas encontradas: %s\n', strjoin(ids, ', '));

    figure('Color','w'); hold on;
    cores = lines(numel(ids));
    for i = 1:numel(ids)
        dados = readmatrix(arquivos{i});
        fprintf('Ks_%s: %d linhas x %d colunas, %d valores NaN\n', ...
                ids{i}, size(dados,1), size(dados,2), nnz(isnan(dados)));

        % pontos originais (plotados mesmo que o ajuste falhe)
        if size(dados,2) >= 2
            plot(dados(:,1), dados(:,2), 'o', 'Color', cores(i,:), ...
                 'MarkerSize', 5, 'HandleVisibility', 'off');
        end

        % curva ajustada (se falhar, avisa qual curva deu problema)
        try
            [~, ~, x_ajuste, y_ajuste] = gerar_ajuste(dados);
            plot(x_ajuste, y_ajuste, '-', 'Color', cores(i,:), 'LineWidth', 2, ...
                 'DisplayName', ['Ajuste: Ks\_' ids{i}]);
        catch ME
            warning('Ajuste de Ks_%s falhou: %s', ids{i}, ME.message);
        end
    end
    grid on; grid minor; box on;
    xlabel('\alpha = b/c');  ylabel('K_s');
    legend('Location','best');
    title('Comparação de Curvas Ajustadas');
    set(gca, 'FontSize', 11, 'LineWidth', 1.1);
    drawnow;
end