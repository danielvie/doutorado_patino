function [y,t,u] = sim_1(config)

    % lendo configuracoes
    Cc    = config.Cc;
    Dc    = config.Dc;
    tstep = config.tstep;
    x0    = config.x0;
    Ts    = config.Ts;

    % simulando com parametros flexiveis
    n_modes = numel(config.modes);
    
    % inicializando saida    
    t = [];
    u = [];
    y = [];
    
    xi0 = x0;
    for i = 1:n_modes
        % lendo modo de operacao (indice do modo inicia em `0`)
        imode = config.modes(i) + 1;
        
        % lendo matrizes A e B 
        Ai = config.Ac{imode};
        Bi = config.Bc{imode};
        
        % calculando ciclo
        ti = (Ts(i):tstep:Ts(i+1)-tstep)';
        ui = ones(size(ti))*config.ur(imode);
        
        yi  = lsim(Ai,Bi,Cc,Dc,ui,ti-ti(1),xi0);
        xi0 = yi(end, :);
        
        % concatenando resultado da iteraco
        t = [t; ti];
        u = [u; ui];
        y = [y; yi];
    end
    
    bla = 1;
end