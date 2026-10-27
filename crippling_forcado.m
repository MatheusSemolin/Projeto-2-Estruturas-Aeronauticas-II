function [sec] = crippling_forcado(sec)

%% Notas:
% Precisa verificar a origem do Sigma_S (não entendi qual valor ele representa nas notas de aula)
% Os valored de F_Ccc e F_ccf provavelmente estão corretos
% Sigma_c também precisa ser conferido, mas deve estar certo
% Subtituir o valor para o K quando o Léo disponibilizar

% Conversões
KSI2MPA = 6.89476;

% Propriedades do material
epslon_0 = 0.002;
FC_y = sec.propriedades.sigma_y;
E = sec.propriedades.E;
K = 0.675;                          % <---- Trocar pelo valor do Léo
t_sk = 0.80;                        % <---- em mm

% Coeficiente C (em ksi !!!)
C = (0.01 * FC_y / KSI2MPA) / (5.88 * sqrt(FC_y / E) + epslon_0);

% Número de painéis
n_panels = size(sec.malha.panels, 1);

% Inicializar vetor de resultados
F_ccf_ksi = zeros(n_panels, 1);
F_ccf_MPa = zeros(n_panels, 1);
R_CF = zeros(n_panels, 1);




%% Parte 1. R_CF - Resistência admissível ao crippling forçado
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_panels

    % Espessura do longeron associada ao painel i
    t_st = sec.malha.panels(i, 3);                          % <---- em mm

    % Resistência admissível ao crippling forçado
    F_ccf_ksi(i) = 20.7 * C * K^(2/3) * (t_st/t_sk)^(1/3);  % <---- em ksi

    F_ccf_MPa(i) = F_ccf_ksi(i) * KSI2MPA;                  % <---- em MPA


end

%% Verificação
if numel(sec.crippling.Fcc_elem) < n_panels
    error(['sec.crippling.Fcc_elem possui %d elemento(s), mas existem ' ...
        '%d painéis em sec.malha.panels.'], ...
        numel(sec.crippling.Fcc_elem), n_panels);
end


%% Razão de crippling forçado:

for i = 1:n_panels

    Sigma_sr_MPa = sec.resultados.sigma_s(i);       % <---- verificar se a origem dos dados está correta
                                                    %       e se as  tensões estão em MPa
    R_CF(i) = Sigma_sr_MPa / F_ccf_MPa(i);          % <---- em MPA

end

% Armazenar resultados
sec.crippling_forcado.F_ccf_ksi = F_ccf_ksi;
sec.crippling_forcado.F_ccf_MPa = F_ccf_MPa;
sec.crippling_forcado.R_CF = R_CF;


%% 2. R_CN - Resistência admissível à compressão normal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Número de painéis
n_panels = size(sec.malha.panels, 1);

% Verificar se existe uma resistência para cada painel
if numel(sec.crippling.Fcc_elem) < n_panels
    error(['sec.crippling.Fcc_elem possui %d elemento(s), mas existem ' ...
        '%d painéis em sec.malha.panels.'], ...
        numel(sec.crippling.Fcc_elem), n_panels);
end

% Verificar se existem tensões suficientes
if numel(sec.resultados.Sigma_nodes) < n_panels
    error(['sec.resultados.Sigma_nodes possui %d elemento(s), mas ' ...
        'existem %d painéis em sec.malha.panels.'], ...
        numel(sec.resultados.Sigma_nodes), n_panels);
end

% Inicializar vetores
F_Ccc_MPa = zeros(n_panels, 1);
R_CN = zeros(n_panels, 1);

for i = 1:n_panels

    % Resistência admissível à compressão normal
    F_Ccc_MPa(i) = sec.crippling.Fcc_elem(i);          % <---- em MPa

    % Tensão normal solicitante [MPa]
    Sigma_c_MPa = sec.resultados.Sigma_nodes(i);       % <---- em MPa

    % Razão entre tensão solicitante e resistência admissível
    R_CN(i) = abs(Sigma_c_MPa) / F_Ccc_MPa(i);                              % <---- confirmar se: 
                                                                            % Sigma_c  é   sec.resultados.Sigma_nodes
                                                                            % F_Ccc    é   ssec.crippling.Fcc_elem

end

% Armazenar resultados
sec.crippling_forcado.F_Ccc_MPa = F_Ccc_MPa;
sec.crippling_forcado.R_CN = R_CN;


%% 3. Cálculo da Margem de Segurança
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inicializar vetores
FS_CF = zeros(n_panels, 1);
MS_CF = zeros(n_panels, 1);

for i = 1:n_panels

    % Combinação das razões de utilização
    FS_CF(i) = 1 / (R_CN(i) + R_CF(i));

    % Margem de segurança
    MS_CF(i) = FS_CF(i) - 1;

end

%% 4. Display dos resultados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:n_panels

    fprintf('=================== Painel %d ===================\n', i);
    fprintf('Sigma_c = %.3f MPa | Sigma_sr = %.3f MPa\n', Sigma_c_MPa, Sigma_sr_MPa);
    fprintf('F_Ccc   = %.3f MPa | F_ccf    = %.3f MPa\n', F_Ccc_MPa(i),  F_ccf_MPa(i));
    fprintf('R_CN    = %.4f      | R_CF     = %.4f\n', R_CN(i), R_CF(i));
    fprintf('FS_CF   = %.4f      | MS_CF    = %.4f\n\n', FS_CF(i), MS_CF(i));

end    

end

