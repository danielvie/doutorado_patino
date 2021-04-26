% [y,t,u] = sim_n(config, Ts)
function [y,t,u] = sim_n(config, Ts, nsim)
    
    config_ = config;
    y   = [];
    t   = [];
    u   = [];
    t0  = 0.0;
    for i = 1:nsim
        % calculo `ek`
        % ek  = x0 - x_target;

        % calculo comando `dtk`
        % dtk = mpc_dualmode_switching(ek,H,Hf,Phi1Np,Qbar,Rbar,Lbar,cbar,Pf,Sf,bf,PhiNp,p);

        % simulando dinamica
        [y_,t_,u_] = sim_1(config_, Ts);
        config_.x0 = y_(end,:)';

        y   = [y;y_];
        t   = [t;t_ + t0];
        u   = [u;u_];
        
        t0  = t(end);
    end
end