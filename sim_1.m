function [y,t,u] = sim_1(config, Ts)
    
    A1 = config.A1;
    A2 = config.A2;
    B1 = config.B1;
    B2 = config.B2;
    Cc = config.Cc;
    Dc = config.Dc;
    
    tstep = config.tstep;
    x0    = config.x0;
    
    ur  = [1, 0];

    % ciclo 1
    t1 = Ts(1):tstep:Ts(2)-tstep;
    u1 = ones(size(t1))*ur(1);
    y1 = lsim(A1,B1,Cc,Dc,u1,t1-t1(1),x0);

    % ciclo 2
    t2 = Ts(2):tstep:Ts(3)-tstep;
    u2 = ones(size(t2))*ur(2);
    y2 = lsim(A2,B2,Cc,Dc,u2,t2-t2(1),y1(end, :));
    
    t  = [t1,t2]';
    u  = [u1,u2]';
    y  = [y1;y2];    
    
end