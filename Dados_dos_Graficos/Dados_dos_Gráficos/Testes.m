clc;

pasta_script = fileparts(mfilename('fullpath'));

cc = 3;
Razao_bc = 0;
Kb = Kb_value(cc, Razao_bc, pasta_script);
fprintf('Valor de Kb para a curva %d na razão %.2f: %.4f\n', cc, Razao_bc, Kb);


arquivos = {fullfile(pasta_script, 'Kb', 'Kb_1.txt'), ...
            fullfile(pasta_script, 'Kb', 'Kb_2.txt'), ...
            fullfile(pasta_script, 'Kb', 'Kb_3.txt'), ...
            fullfile(pasta_script, 'Kb', 'Kb_4.txt')};
plotar_todas_curvas(arquivos);


function [Kb] = Kb_value(cc, Razao_bc, pasta_script)
    % Mapeia a variável cc para o ficheiro txt apontando direto para 'Kb'
    nome_arquivo = fullfile(pasta_script, 'Kb', sprintf('Kb_%d.txt', cc));
    
    if ~isfile(nome_arquivo)
        error('Ficheiro não encontrado: %s\nVerifique se não gravou acidentalmente como "Kb_%d.txt.txt".', nome_arquivo, cc);
    end
    
    dados = readmatrix(nome_arquivo);
    [ftemp] = gerar_ajuste(dados);
    Kb = ftemp(Razao_bc);
end


function [f, ajuste, x_ajuste, y_ajuste] = gerar_ajuste(dados)
    if size(dados,2) ~= 2
        error('A matriz dados deve possuir exatamente duas colunas: [x y].');
    end
    if size(dados,1) < 4
        error('São necessários pelo menos 4 pontos para realizar o ajuste.');
    end
    
    x = dados(:,1);
    y = dados(:,2);
    
    validos = isfinite(x) & isfinite(y);
    x = x(validos);
    y = y(validos);
    
    [x,ordem] = sort(x);
    y = y(ordem);
    
    p = 0.99;
    ajuste = fit(x, y, 'smoothingspline', 'SmoothingParam', p);
    
    f = @(x_consulta) feval(ajuste, x_consulta);
    
    x_ajuste = linspace(min(x),max(x),500);
    y_ajuste = f(x_ajuste);
end


function plotar_todas_curvas(lista_arquivos)
    figure('Color','w');
    hold on;
    
    cores = lines(length(lista_arquivos)); 
    
    for i = 1:length(lista_arquivos)
        caminho_completo = lista_arquivos{i};
        
        if ~isfile(caminho_completo)
            warning('Ficheiro "%s" ignorado pois não foi encontrado.', caminho_completo);
            continue;
        end
        
        dados = readmatrix(caminho_completo);
        x = dados(:,1);
        y = dados(:,2);
        
        [~, ~, x_ajuste, y_ajuste] = gerar_ajuste(dados);
        
        cor = cores(i, :);
        
        [~, nome, ext] = fileparts(caminho_completo);
        nome_legenda = strrep([nome, ext], '_', '\_');
        
        plot(x, y, 'o', 'Color', cor, 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(x_ajuste, y_ajuste, '-', 'Color', cor, 'LineWidth', 2, 'DisplayName', ['Ajuste: ', nome_legenda]);
    end
    
    grid on; 
    grid minor; 
    box on;
    ylabel('K_b');
    xlabel('\alpha = b/c');
    legend('Location','best');
    title('Comparação de Curvas Ajustadas');
    set(gca, 'FontSize', 11, 'LineWidth', 1.1);
    
    hold off;
end