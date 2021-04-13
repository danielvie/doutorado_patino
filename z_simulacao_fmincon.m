%%
% --------------------------------------
% CONFIG INICIAL
% --------------------------------------
clear all
config = init_sim();

%%
% --------------------------------------
% PARAMETROS SIMULACAO
% --------------------------------------

Ts  = [0., 0.2514520, 0.5014520];
dx1 = 0.0;
dx2 = 0.0;

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
x0  = [0.8, 0.5];
lb  = [0.25, 0.25];
ub  = [1.5, 1.5];
A   = [];
b   = [];
Aeq = [];
beq = [];
nonlincon = [];

[x, fval] = fmincon(@(x) fun_custo_patino(config, x), x0, A, b, Aeq, beq, lb, ub, nonlincon, options);

% grafico da trajetoria final
T  = [0, x(1), sum(x)];
yf = sim_1(config, T);
plot(yf(:,1), yf(:,2), 'k', 'linew', 2);

% pontos de `xref` e `x0`
AZUL     = [0.00,0.45,0.74];
VERMELHO = [1.00,0.00,0.00];

plot(config.x0(1), config.x0(2), 'o', 'linew', 2, 'color', AZUL);
h = plot(config.xref(1), config.xref(2), '+', 'linew', 2, 'markersize', 10, 'color', VERMELHO);

grid on;

%
% --------------------------------------
% RESULTADOS
% --------------------------------------

% otimizacao
disp('-----------------------------------')
fprintf('RESULTADOS:\n');
disp(['x    : [', num2str(x, '%.7f, ') ']']);
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

config_ = config;
config_.x0 = config_.x0;
% config_.x0 = [1.873510, -1.118048];

y = sim_n(config_, config_.Ts, nsim);

figure(2);
clf();
hold on;

plot(y(:,1), y(:,2));

plot(config.xref(1), config.xref(2), '+', 'linew', 2, 'markersize', 10, 'color', VERMELHO);


% ultimo ciclo
config_ = config;
config_.x0 = y(end, :);
yend = sim_1(config_, config_.Ts);
plot(yend(:,1), yend(:,2), 'k', 'linew', 2);

fprintf('xend ultimo ciclo: %.6f, %.6f\n', config_.x0(1), config_.x0(2));

%%
% --------------------------------------
% FUNCOES
% --------------------------------------

function Ts = get_ts(x)
    Ts = [0., x(1), sum(x)];
end


function J = fun_custo(config, dT)
    
    T = [0, dT(1), sum(dT)];
    y = sim_1(config, T);

    e = config.x0 - y(end,:);
    J = sum(e.^2);
    
    figure(1);
    hold on;
    plot(y(:,1), y(:,2));
    
end

function J = fun_custo_patino(config, dT)
    
    % lendo parametros de configuracao
    xref = config.xref;
    Q    = config.Q;
    
    % montando vetor de tempo e rodando simulacao
    T  = [0, dT(1), sum(dT)];
    y  = sim_1(config, T);
    
    % grafico trajetorias
    figure(1);
    hold on;
    plot(y(:,1), y(:,2));
        
    % calculando integral erro trajetoria (trapezios)
    ny  = size(y,1);
    I   = 0;
    dx  = config.tstep;
    v_  = 0;
    for i = 1:ny
        % gerando valor da funcao
        xi = (y(i) - xref);
        v  = xi'*Q*xi;

        % atualizando integrador
        I  = I + dx*(v + v_)/2;
        v_ = v;
    end
    
    % calculando distancia x0 xfim
    e  = config.x0 - y(end,:);
    
    % calculando funcao custo
    J  = I + sum(e.^2);
    
end