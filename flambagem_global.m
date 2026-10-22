function [sec] = flambagem_global(sec)
% flambagem_global: Calcula a carga crítica para fambagem global elástica, plástica e induzida por crippling.
%
%   [sec] = flambagem_global(sec)
%
%   ENTRADAS (campos da struct 'sec'):
%       sec.nome      - Identificação da secção (String)
%       sec.E         - Módulo de Elasticidade [MPa]
%       sec.I_eta     - Menor momento de inércia principal [mm^4]
%       sec.L_viga    - Comprimento destravado da coluna [mm]
%       sec.c         - Fator de comprimento efetivo (c) [Adimensional]
%       sec.A         - Área da seccerrtoção transversal [mm^2]
%       sec.N         - Carga axial atuante [N]
%       sec.sigma_E   - Tensão limite elástica do material [MPa]
%
%   SAÍDAS:
%       sec.Cargas_Criticas.Global_elastica - Carga crítica de Euler [N]

E = sec.propriedades.E;
I_eta = sec.propriedades.I_eta;
L_viga = sec.propriedades.L_viga;
c = sec.propriedades.c;
A = sec.propriedades.A;
N = sec.cargas.N;
n = sec.propriedades.n;         
F_07 = sec.propriedades.F_07;  
sigma_y = sec.propriedades.sigma_y;
offset = 0.002;
%%%%%%%%%%%%%%% Flambagem Global Elástica %%%%%%%%%%%%%%
Pcr = (pi^2 * E * I_eta) / (c * L_viga)^2;
if  Pcr/A < sec.propriedades.sigma_E
fprintf('Flambagem Global Elástica (%s) \n', sec.nome)

sec.Cargas_Criticas.Global_elastica = Pcr;
MS_global_elastica = (Pcr/N) - 1;
sec.MS.Global_elastica = MS_global_elastica;

fprintf('   - Carga Crítica (Pcr)     : %.2f N\n', Pcr);
fprintf('   - Margem de Segurança (MS): %.2f\n', MS_global_elastica);
end
%%%%%%%%%%%%%%% Flambagem Global Plástica %%%%%%%%%%%%%%

%%%%%%%%%% Ramberg Osgood usando o método do offset %%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%% Flambagem Global Inelástica %%%%%%%%%%%%%%
if  Pcr/A > sec.propriedades.sigma_E

n = log((3*F_07)/(7*E*offset)) / ...
    log(F_07/sigma_y);    
% Modulo tangente pela equacao de offset
Et = @(sigma) 1 ./ ( ...
    1/E + ...
    (offset*n/sigma_y)*(sigma/sigma_y).^(n-1));

% Equacao de flambagem inelastica
f = @(sigma) sigma*A - ...
    pi^2*I_eta*Et(sigma)/(c*L_viga)^2;

% Tensão crítica elástica de Euler
sigma_euler = Pcr/A;

% Solucao numerica
sigma_cr = fzero(f, [eps, sigma_euler]);

% Carga critica inelastica
Pcr_plastica = sigma_cr*A;

% Margem de seguranca
MS_global_plastica = Pcr_plastica/N - 1;

% Armazenamento
sec.Cargas_Criticas.Global_plastica = Pcr_plastica;

sec.MS.Global_plastica = MS_global_plastica;

% Resultados
fprintf('\nFlambagem Global Inelastica (%s)\n', sec.nome);
fprintf('Expoente n        : %.4f\n', n);
fprintf('Modulo elastico   : %.2f MPa\n', E);
fprintf('Modulo tangente   : %.2f MPa\n', Et(sigma_cr));
fprintf('Tensao critica    : %.2f MPa\n', sigma_cr);
fprintf('Carga critica     : %.2f N\n', Pcr_plastica);
fprintf('Margem seguranca  : %.4f\n', MS_global_plastica);
end


%%%%%%%%%%%%%%% Flambagem Induzida por Crippling %%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Largura Minima Efetiva das Flanges %%%%%%%%%%%%%%%

nodes = sec.malha.nodes;
panels = sec.malha.panels;

% Identificacao das flanges pelos nos das extremidades livres
nos_totais = [panels(:,1);panels(:,2)];
[frequencia, nos_distintos] = groupcounts(nos_totais);
nos_livres = nos_distintos(frequencia == 1);

% Inicializacao
w_min = zeros(length(nos_livres), 1);

for i = 1:length(nos_livres)

    no = nos_livres(i);

    % Encontra o indice do no
    idx = find(nodes(:,1) == no, 1);

    % Tensao de flexao na extremidade livre [MPa]
    fb = sec.resultados.Sigma_nodes(idx);

    % Encontra o painel conectado ao no livre
    p = find(panels(:,1) == no | panels(:,2) == no, 1);

    % Espessura da flange [mm]
    t = panels(p,3);
    if sec.MS.local(p) <0
    
        if fb < 0

            % Modulo da tensao compressiva
            f = abs(fb);

            % Largura minima efetiva [mm]
            w_min(i) = 0.64 * sqrt(E/(2*f)) * t;

            fprintf('\nFlange com extremidade no no %d\n', no);
            fprintf('Tensao de flexao: %.2f MPa\n', fb);
            fprintf('w_min: %.2f mm\n', w_min(i));

        else

            fprintf('\nFlange com extremidade no no %d: tracionada\n', no);

     
        end

    end
end

% Armazena os resultados
sec.crippling.nos_livres = nos_livres;
sec.crippling.w_min = w_min;

Fcc = sec.crippling.Fcc_perfil;
I_eta_crippling = sec.propriedades.I_eta_crippling;
A_crippling = sec.propriedades.A_crippling;
Pcr_crippling_induzido = Fcc * (1 - (Fcc*(c*L_viga)^2)/(pi^2*E*I_eta_crippling)) * A_crippling;
MS_crippling_induzido = Pcr_crippling_induzido/N - 1;
sec.MS.crippling_induzido = MS_crippling_induzido;


sec.Cargas_Criticas.Pcr_crippling_induzido = Pcr_crippling_induzido;





end