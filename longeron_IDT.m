
function [sec] = longeron_IDT(sec)

%% Cálculo da tensão semi-diagonal do longeron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KSI2MPA = 6.89476;


% Coefientes da IDT
alpha = 45;                          % <---- Trocar pelo valor do Léo
k = 0.675;                              % <---- Trocar pelo valor do Léo

% Número de painéis
n_panels = size(sec.malha.panels, 1);

% Inicializar vetores
A = zeros(n_panels, 1);
B = zeros(n_panels, 1);
R_T = zeros(n_panels, 1);
F_IDT = zeros(n_panels, 1);

for i = 1:n_panels

% Entradas necessárias

    f_sk = NaN;       % Tensão cisalhante no painel
    f_cl = NaN;       % Compressão no longeron
    f_ccr_0 = NaN;    % Resistência crítica em compressão
    f_scr_0 = NaN;    % Resistência crítica em cisalhamento


% Grandezas geométricas necessárias]

    A_st = NaN;       % Área da seção do stiffener
    h = NaN;          % Dimensão geométrica conforme a referência
    t = NaN;          % Espessura conforme a referência


% Cálculo


    A(i) = f_scr_0 / f_ccr_0;
    B(i) = f_sk / f_cl;
    
    R_T(i) = -A(i)/(2*B(i)) + ...
             sqrt((A(i)/(2*B(i)))^2 + 1);
    
    F_IDT(i) = (k*f_sk*cot(alpha)) / ...
               (A_st/(h*t) + 0.5*(1-k)*R_T(i));

end

% Armazenar resultados
sec.longeron_IDT.A = A;
sec.longeron_IDT.B = B;
sec.longeron_IDT.R_T = R_T;
sec.longeron_IDT.F_IDT = F_IDT;

end
