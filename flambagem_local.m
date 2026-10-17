function [sec] = flambagem_local(sec)
% Constantes
E = sec.propriedades.E; % [MPa] - Incluí o Módulo de Elasticidade na fórmula!
v = sec.propriedades.v; % Coeficiente de Poisson
N = sec.cargas.N;
A_total = sec.propriedades.A;
nodes  = sec.malha.nodes;
panels = sec.malha.panels;

b = zeros(size(panels,1), 1);
sigma_b = zeros(size(panels,1),1);
sigma_s = zeros(size(panels,1),1);
for p = 1:size(panels,1)
    n1 = panels(p,1);  
    n2 = panels(p,2); 
    % largura do painel 
    xy1 = nodes(nodes(:,1) == n1, 2:3);
    xy2 = nodes(nodes(:,1) == n2, 2:3);
    b(p) = hypot(xy2(1) - xy1(1), xy2(2) - xy1(2));

    idx1 = find(nodes(:,1) == n1);
    idx2 = find(nodes(:,1) == n2);

    % Puxa as tensões normais calculadas na função Normal(sec)
    sig1 = sec.resultados.Sigma_nodes(idx1);
    sig2 = sec.resultados.Sigma_nodes(idx2);

    % Guarda o valor mais negativo (maior compressão)
    sigma_b(p) = min(sig1, sig2);
    
    if sigma_b(p) >=0
        sigma_b(p) = NaN;
    end

    tau_tor = sec.resultados.tau_torcao(p); 

    % Puxa a função do tau fletor específica deste painel
    fun_tau_flex = sec.resultados.Tau_fletor_sym{p}; 

    % Integra numericamente de 0 a L e divide por L
    tau_flex_medio = integral(fun_tau_flex, 0, b(p)) / b(p);

    % Soma os efeitos (atenção para garantir que os sinais/sentidos estejam coerentes)
    sigma_s(p) = abs(tau_tor + tau_flex_medio);
end

% Extrai os dados geométricos como vetores coluna (N x 1)
t = sec.malha.panels(:,3);

% Monta a matriz K (N x 3) unindo os vetores
% Se já estiver tudo numa matriz, basta usar K = sec.sua_matriz_K;
K_matriz = [sec.coef_K.Kc, sec.coef_K.Kb, sec.coef_K.Ks];

% [sigma_cr_c, sigma_cr_b, sigma_cr_s]
sigma_cr_matriz = (K_matriz .* (t ./ b).^2 * pi^2 * E) / (12 * (1 - v^2));

R_c = abs((N/A_total)./sigma_cr_matriz(:,1));
R_b = abs(sigma_b./sigma_cr_matriz(:,2));
R_s = abs(sigma_s./sigma_cr_matriz(:,3));


% Margens de Segurança MS (perguntar pro Murilo)
% mask_b = (R_b > 0.1);
% mask_c = (R_c > 0.1);
% mask_s = (R_s > 0.1);

 mask_b = ~isnan(R_b);
 mask_c = ~isnan(R_c);
 mask_s = ~isnan(R_s);

caso_bc = mask_b & mask_c & ~mask_s;
caso_cs = ~mask_b & mask_c & mask_s;
caso_bcs = mask_b & mask_c & mask_s;
caso_bs = mask_b & ~mask_c & mask_s;
% Perguntar pro Murilo
caso_b = mask_b & ~mask_c & ~mask_s;
caso_c = ~mask_b & mask_c & ~mask_s;
caso_s = ~mask_b & ~mask_c & mask_s;

% 3. Inicializa o vetor de Margem de Segurança (MS) preenchido com NaN
MS = nan(size(R_b));
MS(caso_bc)  = MS_bc(R_b(caso_bc), R_c(caso_bc));      
MS(caso_cs)  = MS_cs(R_c(caso_cs), R_s(caso_cs));
MS(caso_bs)  = 1 ./ sqrt(R_b(caso_bs).^2 + R_s(caso_bs).^2) - 1;
MS(caso_bcs)  = 1 ./ sqrt(R_b(caso_bcs).^2 +  R_c(caso_bcs).^2 + R_s(caso_bcs).^2) - 1;
% Perguntar pro Murilo
MS(caso_b) = 1 ./ R_b(caso_b) - 1;
MS(caso_c) = 1 ./ R_c(caso_c) - 1;
MS(caso_s) = 1 ./ R_s(caso_s) - 1;

sec.MS.local = MS;
sec.R.R_c = R_c;
sec.R.R_b = R_b;
sec.R.R_s = R_s;


sec.resultados.sigma_cr_Kc = sigma_cr_matriz(:, 1);
sec.resultados.sigma_cr_Kb = sigma_cr_matriz(:, 2);
sec.resultados.sigma_cr_Ks = sigma_cr_matriz(:, 3);
end

function MS = MS_bc(Rb, Rc)
% Caso flexão + compressão: Rc + Rb^1.75 = 1.
% Em B: lambda*Rc + (lambda*Rb)^1.75 = 1 (sem solução fechada, usa fzero).
% Os dois termos crescem com lambda, então lambda <= 1/max(Rb,Rc).
MS = nan(size(Rb));
for i = 1:numel(Rb)
    g = @(lam) lam*Rc(i) + (lam*Rb(i)).^1.75 - 1;
    lam = fzero(g, [0, 1/max(Rb(i), Rc(i))]);
    MS(i) = lam - 1;
end
end

function MS = MS_cs(Rc, Rs)
% Caso compressão + cisalhamento: Rc + Rs^2 = 1.
% Em B: lambda*Rc + lambda^2*Rs^2 = 1  ->  lambda = 2/(Rc + sqrt(Rc^2 + 4*Rs^2))
lambda = 2 ./ (Rc + sqrt(Rc.^2 + 4*Rs.^2));
MS = lambda - 1;
end

