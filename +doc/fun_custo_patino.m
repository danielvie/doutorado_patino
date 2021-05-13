function J = fun_custo_patino(config, X)
    
    % lendo parametros de configuracao
    xref = config.xref;
    Q    = config.Q;
    
    nmodes = numel(config.modes);
    nstate = numel(config.xref);
    
    dT   = X(1:nmodes);
    x0   = X(nmodes+1:nmodes+nstate);
    
    % ajustando condicao inicial
    config_ = config;
    config_.x0 = x0;
    
    % montando vetor de tempo e rodando simulacao
    Ts = zeros(1, numel(dT) + 1);
    for i = 2:numel(Ts)
        Ts(i) = Ts(i-1) + dT(i-1);
    end
    
    config_.Ts = Ts;
    y  = doc.sim_1(config_);
    
    y_ = y - y(1,:);
    
    if size(y_, 2) > 2
        figure(4);
        plot(y_(:,1), y_(:,2)); hold on;
        plot(y_(1,1), y_(1,2), 's');
        plot(y_(end,1), y_(end,2), 'o');

        axis([-0.5, 1.5, -0.5, 1]);
        title('Tensao C1 vs Tensao C2')
        grid on;
        hold off;

        figure(5);
        plot(y_(:,1), y_(:,3)); hold on;
        plot(y_(1,1), y_(1,3), 's');
        plot(y_(end,1), y_(end,3), 'o');

        axis([-0.5, 1.5, -0.5, 0.5]);
        title('Tensao C1 vs Corrente Il')
        grid on;
        hold off;

        figure(6);
        plot(y_(:,2), y_(:,3)); hold on;
        plot(y_(1,2), y_(1,3), 's');
        plot(y_(end,2), y_(end,3), 'o');

        axis([-0.5, 1.0, -0.3, 0.3]);
        title('Tensao C2 vs Corrente Il')
        grid on;
        hold off;
    end
    
    % axis([-0.1, 2, -0.1, 2, -0.2, 0.2]);
    
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
    
    % saida resultado funcao
    J = I;    
end