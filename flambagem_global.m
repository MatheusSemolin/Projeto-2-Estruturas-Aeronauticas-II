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
%%%%%%%%%%%%%%% Flambagem Global Elástica %%%%%%%%%%%%%%
Pcr = (pi^2 * E * I_eta) / (c * L_viga)^2;
fprintf('1. Flambagem Global Elástica (%s) \n', sec.nome)
if  Pcr/A > sec.propriedades.sigma_E
    fprintf("Flambagem global Elástica não ocorre\n")
end
sec.Cargas_Criticas.Global_elastica = Pcr;
MS_global_elastica = (Pcr/N) - 1;
sec.MS.Global_elastica = MS_global_elastica;

fprintf('   - Carga Crítica (Pcr)     : %.2f N\n', Pcr);
fprintf('   - Margem de Segurança (MS): %.2f\n', MS_global_elastica);
%%%%%%%%%%%%%%% Flambagem Global Plástica %%%%%%%%%%%%%%
syms sigma_sym

k = (3/7) *(E/F_07)^(n-1); % Conversão realizada com base na equação de Ramberg-Osgood do livro do Bruhn e das fórmulas anotadas em sala.

E_t = ((1/E)+ (n*k*sigma_sym^(n-1))/E^n)^-1;
f_sym = sigma_sym*A - (E_t*I_eta*pi^2)/(c*L_viga)^2;
df_sym = diff(f_sym,sigma_sym);

f_num = matlabFunction(f_sym);
df_num = matlabFunction(df_sym);

sigma_cr = sec.propriedades.sigma_E;
for i = 1:100
    
    if abs(f_num(sigma_cr)) < 1e-6
        break;
    end
    sigma_cr = sigma_cr - f_num(sigma_cr)/df_num(sigma_cr);
end

Pcr_plastica = sigma_cr * A;
MS_global_plastica = (Pcr_plastica / N) - 1;
sec.MS.Global_plastica = MS_global_plastica;
sec.Cargas_Criticas.Global_plastica = Pcr_plastica;
fprintf('\n2. Flambagem Global Plástica (%s) \n', sec.nome);
    fprintf('   - Carga Crítica (Pcr)     : %.2f N\n', Pcr_plastica);
    fprintf('   - Margem de Segurança (MS): %.2f\n', MS_global_plastica);
    fprintf('============================================================\n\n');
end