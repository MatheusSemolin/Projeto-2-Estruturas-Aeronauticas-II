clc
clear variables

% =========================================================================
% MARGENS DE SEGURANCA DOS LONGERONS A CRIPPLING FORCADO
% =========================================================================
clear; clc;

KSI2MPA = 6.89476;
tau_apl = 72.000;
h = 185.000;
t = 0.800;

alpha_deg = 46.557;
alpha = deg2rad(alpha_deg);
K = 0.465;

% 2. DADOS DOS LONGERONS (Seccao J e Seccao Z)
nomes = ["Longeron J", "Longeron Z"];
Area_real = [233.344, 135.708];
N_cargas = [-1560.000, -1560.000];
t_bar = [1.172, 1.000];
F_cc_cr_MPa = [194.272, 233.793];

% 3. PROPRIEDADES DO MATERIAL (Al 7075-T6)
sigma_y = 468.843;
E_st = 72394.980;
eps_0 = 0.002;

Fcy_ksi = sigma_y / KSI2MPA;
C_mat = (0.01 * Fcy_ksi) / (5.88 * sqrt((sigma_y / E_st) + eps_0));

fprintf('======================================================\n');
fprintf(' RESULTADOS - CRIPPLING FORCADO (LONGERONS)\n');
fprintf('======================================================\n\n');

for i = 1:2
    nome = nomes(i);
    A = Area_real(i);
    N = abs(N_cargas(i));
    
    sigma_induzida = tau_apl * ((h * t) / A) * cot(alpha);
    sigma_c_aplicada = N / A;
    
    F_cc_f_ksi = 20.7 * C_mat * (K^(2/3)) * ((t_bar(i) / t)^(1/3));
    F_cc_f_MPa = F_cc_f_ksi * KSI2MPA;
    
    R_cn = sigma_c_aplicada / F_cc_cr_MPa(i);
    R_cf = sigma_induzida / F_cc_f_MPa;
    
    MS = (1 / (R_cn + R_cf)) - 1;
    
    fprintf('--- %s ---\n', nome);
    fprintf('Tensao Compressao: %8.2f MPa\n', sigma_c_aplicada);
    fprintf('Tensao Induzida  : %8.2f MPa\n', sigma_induzida);
    fprintf('F_cc,cr          : %8.2f MPa\n', F_cc_cr_MPa(i));
    fprintf('F_cc,f           : %8.2f MPa\n', F_cc_f_MPa);
    fprintf('R_cn             : %8.4f\n', R_cn);
    fprintf('R_cf             : %8.4f\n', R_cf);
    fprintf('Margem de Seguranca: %8.4f\n', MS);
    
    if MS > 0
        fprintf('ESTADO: SEGURO\n\n');
    else
        fprintf('ESTADO: FALHA ESTRUTURAL\n\n');
    end
end
