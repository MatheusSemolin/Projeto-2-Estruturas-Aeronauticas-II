%% SCRIPT PRINCIPAL - ANÁLISE DA SEÇÃO C
clc; clear variables; close all;

% Nós: [ID, x, y, Area_do_Boom]
sec.malha.nodes = [
    1,  0,  1.5, 22.811; 
    2, 58.5,  1.5,0;
    3, 58.5, 55, 60;
    4, 58.5, 78.5, 0;
    5, 0, 78.5, 22.811;
    ];

% Espessuras das chapas
t1 = 3; % [mm]

% Painéis: [ID_Nó_Origem, ID_Nó_Destino, Espessura_da_Chapa]
sec.malha.panels = [
    1,  2,  t1; 
    2,  3,  t1; 
    3,  4,  t1; 
    4,  5,  t1; 
    ];

% =========================================================================
% Cargas
% =========================================================================
sec.cargas.Vx = 3.5e3; %N
sec.cargas.Vy = 18e3; %N
sec.cargas.Mx = -37e3; %N.mm
sec.cargas.My = 115e3; %N.mm
sec.cargas.Mz = 130e3; %N.mm
sec.cargas.N  = 37.6e3;  %N.mm     % Carga axial atuante [N]

% =========================================================================
% Propriedades
% =========================================================
sec.nome      = 'Perfil Caverna C';
sec.propriedades.Ixx = 701139.04;
sec.propriedades.Iyy = 305215.75;
sec.propriedades.Ixy = 18706.97;
sec.propriedades.centroide = [40.215, 41.310];
sec.propriedades.centro_de_cisalhamento = [85.75, 41.97];
sec.propriedades.I_eta     = 304333.83;      % Menor momento de inércia principal [mm^4]
sec.propriedades.A         = 687.09;         % Área da secção transversal [mm^2]
sec.propriedades.E         = 72395;       % Módulo de Elasticidade [MPa]
sec.propriedades.v = 0.33; % Coeficiente de Poisson
sec.propriedades.L_viga    = 185;         % Comprimento destravado da coluna [mm]
sec.propriedades.c         = 2.0;         % Fator de comprimento efetivo (c) [Adimensional]
sec.propriedades.sigma_E   = 370;         % Tensão limite elástica do material [MPa]
sec.propriedades.sigma_y   = 468.843;     % Tensão de escoamento à compresão [MPa]
sec.propriedades.n         = 9.2;         % Expoente de Ramberg-Osgood
sec.propriedades.F_07      = 482;         % MPa Tensão de referência [MPa]

sec.propriedades.A_crippling = sec.propriedades.A; %% Não tem flambagem local então não tem w min
sec.propriedades.I_eta_crippling = sec.propriedades.I_eta; %% Não tem flambagem local então não tem w min

fprintf('\n====================================================================\n');
fprintf(' ANÁLISE: %s\n', sec.nome);
fprintf('====================================================================\n');


% =========================================================================
% Calculo de Tensões
% =========================================================================
[sec] = Normal(sec);
[sec] = cisalhamento_aberto(sec);
[sec] = Cisalhamento_torcor(sec);


% =========================================================================
%  CRIPPLING NORMAL
% =========================================================================
sec = crippling_normal(sec);

% =========================================================================
%  FLAMBAGEM LOCAL
% =========================================================================
fprintf('\n---------------------------------------------------------\n');
fprintf(' >> ANÁLISE DE FLAMBAGEM LOCAL\n');
fprintf('---------------------------------------------------------\n');sec = coeficientes_paineis(sec,'apoiado');

sec = flambagem_local(sec);
T = table((1:size(sec.malha.panels,1))', sec.R.R_c, sec.R.R_b, ...
    sec.R.R_s, sec.MS.local, ...
    'VariableNames', {'Painel','R_c','R_b','R_s','MS'});
disp(T);
disp('Processo concluído com sucesso!');

% =========================================================================
%  FLAMBAGEM GLOBAL (Elástica e Plástica)
% =========================================================================
fprintf('\n---------------------------------------------------------\n');
fprintf(' >> ANÁLISE DE FLAMBAGEM GLOBAL \n');
fprintf('---------------------------------------------------------\n');
[sec] = flambagem_global(sec);

disp('_____________________________________________________________');


% =========================================================================
%  EXIBIÇÃO DOS RESULTADOS FINAIS 
% =========================================================================
N_paineis = length(sec.MS.local);

Modo_de_Falha = [
    {'Flambagem Global (Elástica)'};
    {'Flambagem Global (Plástica)'};
    {'Flambagem Global (Crippling Induzido)'}; 
    {'Crippling Normal'}; 
    compose('Flambagem Local (Painel %d)', (1:N_paineis)')
    ];

Margem_de_Seguranca = [
    sec.MS.Global_elastica; 
    sec.MS.Global_plastica;
    sec.MS.crippling_induzido; 
    sec.MS.crippling_normal; 
    sec.MS.local(:) 
    ];

Tabela_MS = table(Modo_de_Falha, Margem_de_Seguranca, 'VariableNames', {'Modo_de_Falha', 'MS'});

fprintf('\n================ RESULTADOS FINAIS ================\n');
fprintf('PARTE A - REFORÇADORES\n');
fprintf('Margens de segurança à flambagem e instabilidades\n\n');
disp(Tabela_MS);
fprintf('===================================================\n');

fprintf("Obs: Os valores que aparecem 'NaN' representam quando a condição não se aplica, não ocorre ou a margem não existe \n")