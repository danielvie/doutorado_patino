%%
% --------------------------------------
% CONFIG INICIAL
% --------------------------------------
clear all
config = init_sim();

cores.AZUL     = [0.00,0.45,0.74];
cores.VERMELHO = [1.00,0.00,0.00];

%%
% --------------------------------------
% PARAMETROS SIMULACAO
% --------------------------------------

Ts  = [0., 0.2514520, 0.5014520];
dx1 = 0.;
dx2 = 0.;

config.Ts    = Ts;
config.x0    = [1.870801, -1.119853] + [dx1,dx2];
% config.x0    = [1.870801, -1.119853];
config.tstep = 1e-5;

disp('-----------------------------------');
fprintf('x0: %.6f, %.6f\n', config.x0(1), config.x0(2));
disp('-----------------------------------');

%%
% --------------------------------------
% FMINCON
% --------------------------------------

figure(1);
clf();

% 'Algorithm': 'active-set', 'interior-point', 'sqp', 'sqp-legacy', 'trust-region-reflective'
% options = optimoptions('fmincon', 'display', 'iter' , 'DiffMinChange', config.tstep*10, 'Algorithm', 'sqp');
options = optimoptions('fmincon', 'DiffMinChange', config.tstep*10, 'Algorithm', 'sqp');
% options = optimoptions('fmincon', 'Algorithm', 'active-set');
x0  = [0.25, 0.25, config.x0];
lb  = [0.25, 0.25, config.x0 - 2.0];
ub  = [1.50, 1.50, config.x0 + 2.0];
A   = [];
b   = [];
Aeq = [];
beq = [];
nonlincon = [];

[x, fval] = fmincon(@(x) fun_custo_patino(config, x), x0, A, b, Aeq, beq, lb, ub, nonlincon, options);
% [x, fval] = fmincon(@(x) fun_custo(config, x), x0, A, b, Aeq, beq, lb, ub, nonlincon, options);

% grafico da trajetoria final
T  = [0, x(1), sum(x(1:2))];

config_    = config;

config_.Ts = get_ts(x(1:2));
config_.x0 = x(3:4);
yf = sim_1(config_, T);
plot(yf(:,1), yf(:,2), 'k', 'linew', 2);

% pontos de `xref` e `x0`

plot(config_.x0(1), config_.x0(2), 'o', 'linew', 2, 'color', cores.AZUL);
h = plot(config_.xref(1), config_.xref(2), '+', 'linew', 2, 'markersize', 10, 'color', cores.VERMELHO);

grid on;

%
% --------------------------------------
% RESULTADOS
% --------------------------------------

% otimizacao
disp('-----------------------------------')
fprintf('RESULTADOS:\n');
disp(['x    : [', num2str(x, '%.7f, ') ']']);
disp(['dT   : [', num2str(x(1:2), '%.7f, ') ']']);
disp(['Ts   : [', num2str(get_ts(x), '%.7f, ') ']']);
disp(['fval : ', num2str(fval)]);

% pontos de cruzamento com as linhas de `xref`
% [Cx, Cy] = cruzamentos(yf, config.xref);
xmean = mean(yf);

diffx = config.xref(1) - xmean(1);
diffy = config.xref(2) - xmean(2);

fprintf('diferença centro cruzamento e `xref` : [%.6f, %.6f]\n', diffx, diffy);
disp(['fval : ', num2str(fval)]);

%%
% --------------------------------------
% SIMULACAO COMPLETA
% -------------------------------------- 
nsim = 50;


disp(config_.x0);
y = sim_n(config_, config_.Ts, nsim);


figure(2);
clf();
hold on;

plot(y(:,1), y(:,2));

plot(config.xref(1), config.xref(2), '+', 'linew', 2, 'markersize', 10, 'color', cores.VERMELHO);

% ultimo ciclo
c    = config_;
c.x0 = y(end,:);
yend = sim_1(c, c.Ts);

% ymean = mean(yend);
% plot(ymean(1), ymean(2), 'o', 'linew', 2, 'markersize', 10, 'color', cores.VERMELHO);

plot(yend(:,1), yend(:,2), 'k', 'linew', 2);
axis equal;
fprintf('xend ultimo ciclo: %.6f, %.6f\n', config_.x0(1), config_.x0(2));
grid on;

%%
% --------------------------------------
% FUNCOES
% --------------------------------------

function Ts = get_ts(x)
    x_ = x(1:2);
    Ts = [0., x_(1), sum(x_)];
end

function J = fun_custo(config, X)
    
    % lendo parametros de configuracao
    dT = X(1:2);
    x0 = X(3:4);
    
    % montando vetor de tempo e rodando simulacao
    config_ = config;
    config_.x0 = x0;
    
    Ts = get_ts(dT);
    y  = sim_1(config_, Ts);
    
    % calculando erro
    e  = config_.x0 - y(end,:);
    
    % calculando erro com xref
    ymean = mean(y);
    eref  = config_.xref - ymean;
    
    % calculando valor de custo
    % TODO: ajustar estimativa de `eref` para usar no custo
    J  = sum(e.^2) + sum(eref.^2)*0.05;
    
    % grafico de trajetoria
    figure(1);
    hold on;
    plot(y(:,1), y(:,2));
    title(sprintf('diff x0: %.4f, %.4f; dT: %.4f, %.4f', e(1), e(2), dT(1), dT(2)));
end

function J = fun_custo_patino(config, X)
    
    % lendo parametros de configuracao
    xref = config.xref;
    Q    = config.Q;
    
    dT   = X(1:2);
    x0   = X(3:4);
    
    % ajustando condicao inicial
    config_ = config;
    config_.x0 = x0;
    
    % montando vetor de tempo e rodando simulacao
    T  = [0, dT(1), sum(dT)];
    y  = sim_1(config_, T);
    
    % grafico trajetorias
    figure(1);
    hold on;
    plot(y(:,1), y(:,2));
    
        
    % calculando integral erro trajetoria (trapezios)
    ny  = size(y,1);
    I   = 0;
    dx  = config_.tstep;
    v_  = 0;
    for i = 1:ny
        % gerando valor da funcao
        xi = (y(i,:) - xref);
        v  = xi*Q*xi';

        % atualizando integrador
        I  = I + dx*(v + v_)/2;
        v_ = v;
    end
    
    % calculando distancia x0 xfim
    e  = config_.x0 - y(end,:);
    
    % calculando erro com xref
    ymean = mean(y);
    eref  = config_.xref - ymean;
    
    % calculando funcao custo
    C  = [I, sum(e.^2), sum(eref.^2)]; % elemantos para soma custo
    Cn = C.*[1, 10000, 100]; % custo normalizado
    J  = sum(Cn); % soma custos
        
    title(sprintf('diff x0: %.4f, %.4f; dT: %.4f, %.4f', e(1), e(2), dT(1), dT(2)));
end