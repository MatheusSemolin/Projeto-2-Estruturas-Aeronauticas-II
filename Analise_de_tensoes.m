%% SCRIPT PRINCIPAL - ANÁLISE COMPLETA DA SEÇÃO (Cisalhamento + Flambagem Global)
clc; clear variables; close all;

% Nós: [ID, x, y, Area_do_Boom]
sec.malha.nodes = [
    1,  0,  10.8, 0; 
    2,  0,  0, 14.811;
    3,   40.60, 0.4, 0;
    4, 40.60, 51.2, 2*14.811;
    5, 10, 51.2, 0;
    6, 72, 51.2, 0
];

% Espessuras das chapas
t1 = 0.8;
t2 = 2;

% Painéis: [ID_Nó_Origem, ID_Nó_Destino, Espessura_da_Chapa]
sec.malha.panels = [
    1,  2,  t1; 
    2,  3,  t1; 
    3,  4,  t2; 
    4,  5,  t1; 
    4,  6,  t1;
];

% =========================================================================
% Cargas
% =========================================================================
sec.cargas.Vx = 1150;
sec.cargas.Vy = 860; 
sec.cargas.Mx = -2400;
sec.cargas.My = 85000;
sec.cargas.Mz = 400;
sec.cargas.N  = 1560;       % Carga axial atuante [N]

% =========================================================================
% Propriedades
% =========================================================
sec.nome      = 'Perfil Longarina J';
sec.propriedades.Ixx = 99347.03;
sec.propriedades.Iyy = 42829.78;
sec.propriedades.Ixy = 26755.18;
sec.propriedades.centroide = [36.694, 28.395];
sec.propriedades.centro_de_cisalhamento = [47.10, 22.59];
sec.propriedades.I_eta     = 32173.24;      % Menor momento de inércia principal [mm^4]
sec.propriedades.A         = 237.42;         % Área da secção transversal [mm^2]
sec.propriedades.E         = 72395;       % Módulo de Elasticidade [MPa]
sec.propriedades.L_viga    = 320;         % Comprimento destravado da coluna [mm]
sec.propriedades.c         = 1.0;         % Fator de comprimento efetivo (c) [Adimensional]
sec.propriedades.sigma_E   = 370;         % Tensão limite elástica do material [MPa]
sec.propriedades.n         = 9.2;         % Expoente de Ramberg-Osgood
sec.propriedades.F_07      = 482;         % MPa Tensão de referência [MPa]

% =========================================================================
% Calculo de Tensões
% =========================================================================
%[sec] = Normal(sec);
[sec] = cisalhamento_aberto(sec);
%[sec] = Cisalhamento_torcor(sec);


% =========================================================================
%  FLAMBAGEM GLOBAL (Elástica e Plástica)
% =========================================================================
disp('Iniciando análise de flambagem global...');
[sec] = flambagem_global(sec);

disp('---------------------------------------------------------');
disp('Processo concluído com sucesso!');

%%


% 2. Junta tudo em um único vetor de nós
todos_os_nos = [sec.malha.panels(:, 1); sec.malha.panels(:, 2)];

% 3. Encontra nós únicos e conta quantas vezes cada um aparece na malha
[nos_unicos, ~, ic] = unique(todos_os_nos);
contagens = accumarray(ic, 1);

% 4. Filtra os nós que aparecem 2 ou mais vezes (Nós compartilhados / não-livres)
nos_compartilhados = nos_unicos(contagens >= 2);

% 5. Filtra os nós que aparecem APENAS 1 vez (Pontas livres / Flanges livres)
nos_livres = nos_unicos(contagens == 1);

% Exibe no Command Window para conferir
disp('Nós compartilhados (Junções / Almas):');
disp(nos_compartilhados');

disp('Nós de extremidade (Pontas livres / Flanges):');
disp(nos_livres');


