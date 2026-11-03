% =========================================================================
% MARGENS DE SEGURANCA A TENSAO SEMI-DIAGONAL
% =========================================================================
clear; clc;

% 1. DADOS DIRETOS DO EXCEL
tau_apl = 72.000;             % Tensao aplicada [MPa]
d = 320.000;                  % Distancia longitudinal [mm]
h = 185.000;                  % Comprimento de arco [mm]
t = 0.800;                    % Espessura da chapa [mm]

% Valores convergidos retirados do Excel
sigma_s_cr_0 = 14.006;        % f_s_cr [MPa]
alpha_deg = 46.529;           % [deg]
alpha = deg2rad(alpha_deg);   % [rad]
K = 0.465;                    % Fator K do Excel 

% 2. DADOS DOS REFORCADORES (Seccao J, Seccao Z, Caverna C)
nomes = ["Longeron J", "Longeron Z", "Caverna C"];
is_longeron = [1, 1, 0]; 

% ÁREAS E CARGAS
Area_real = [233.344, 135.708, 687.073]; % [mm^2]
N_cargas  = [-1560.000, -1560.000, -37600.000]; % Cargas axiais [N]
n_ref = [2.000, 2.000, 2.000];           % Consideracao de simetria (2 por painel)

% TENSÕES CRÍTICAS E ADMISSÍVEIS[cite: 1]
s_c_cr_0 = [11.460, 11.460, 0.000];      % f_c_cr [MPa]
s_adm = [235.427, 153.315, 0.000];       % Tensao admissivel [MPa]

% Parâmetros da Caverna C (Calculados dinamicamente com o K fixo)
F_cc_cr_MPa = 239.459; 
sigma_y = 500.000; eps_y = 0.007; E_st = 72394.980; eps_0 = 0.002;
C_mat = eps_y / (5.88 * sqrt((sigma_y/E_st) + eps_0));
t_bar_cav = 3.000; 
F_cc_f_ksi = 20.7 * C_mat * (K^(2/3)) * (t_bar_cav / t);
F_cc_f_MPa = F_cc_f_ksi * 6.895;         % Tensão de crippling corrigida [MPa]

fprintf('======================================================\n');
fprintf(' RESULTADOS DAS MARGENS DE SEGURANCA\n');
fprintf('======================================================\n\n');

for i = 1:3
    nome = nomes(i); 
    A = Area_real(i);
    N = abs(N_cargas(i)); % Tira o sinal para nao quebrar as contas de flambagem

    if is_longeron(i)
        % Tensao induzida pelo campo de tracao
        sigma_induzida = tau_apl * ((h*t)/A) * (1/tan(alpha)); 

        % Simetria: Carga Dobrada e 2 Reforçadores
        P_axial = 2 * N;
        Area_colaborante = (n_ref(i) * A) + (h * t);
        sigma_c = P_axial / Area_colaborante;

        % Parametros de interaçao
        param_A = sigma_s_cr_0 / s_c_cr_0(i);
        param_B = tau_apl / sigma_c;

        R_tau = -(param_A / (2*param_B)) + sqrt((param_A^2 / (4*param_B^2)) + 1);
        R_c = (param_A / param_B) * R_tau;
        sigma_c_cr = R_c * s_c_cr_0(i);

        D = 0.445 * ((h*t)/(n_ref(i) * A)) * R_c * sqrt(abs(sigma_c_cr)); 
        termo_raiz = max(0, D^2 - sigma_c);

        sigma_c_bar = sigma_c + 2*D*sqrt(termo_raiz);
        sigma_apl = sigma_c_bar + sigma_induzida; 

        MS = (s_adm(i) / sigma_apl) - 1; 
    else
        % Caverna C
        sigma_induzida = tau_apl * ((d*t)/A) * tan(alpha);
        sigma_c_aplicada = (N / A); 

        R_cu = sigma_c_aplicada / F_cc_cr_MPa;
        R_cf = sigma_induzida / F_cc_f_MPa;

        MS = (1 / (R_cu + R_cf)) - 1;
    end

    fprintf('--- %s ---\n', nome);
    fprintf('Tensao Induzida: %8.2f MPa\n', sigma_induzida);
    fprintf('Margem de Seguranca: %8.4f\n', MS);

    if MS > 0
        fprintf('ESTADO: SEGURO\n\n');
    else
        fprintf('ESTADO: FALHA ESTRUTURAL\n\n');
    end
end
