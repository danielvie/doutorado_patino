%%
% --------------------------------------
% CONFIG INICIAL
% --------------------------------------
clear all
config = init_sim_2();

cores.AZUL     = [0.00,0.45,0.74];
cores.VERMELHO = [1.00,0.00,0.00];

%%
% --------------------------------------
% PARAMETROS SIMULACAO
% --------------------------------------


%%
% --------------------------------------
% SIMULACAO COMPLETA
% -------------------------------------- 
nsim = 50;

config_ = config;

config_.x0 = config.xref;
% [y,t,u] = sim_1(config_);

[y, t] = sim_n(config_, nsim);

close all;
plot(t, y(:,3));


%%
% --------------------------------------
% FUNCOES
% --------------------------------------

function Ts = get_ts(x)
    x_ = x(1:2);
    Ts = [0., x_(1), sum(x_)];
end

function [c,ceq] = nonlincon_patino(config, X)
    % c(x)   <= 0
    % xeq(x)  = 0
    
    % lendo parametros de configuracao
    dT   = X(1:2);
    x0   = X(3:4);
    
    % ajustando condicao inicial
    config_ = config;
    config_.x0 = x0;
    
    % montando vetor de tempo e rodando simulacao
    T  = [0, dT(1), sum(dT)];
    y  = sim_1(config_, T);
    
    e  = config_.x0 - y(end,:);
    
    c   = [];
    
    ceq = sum(e.^2);
    
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
    dx  = config_.tstep;
    I   = 0;
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

    J = I;
    title(sprintf('diff x0: %.7f, %.7f; dT: %.7f, %.7f; J: %.7f', e(1), e(2), dT(1), dT(2), J));
end