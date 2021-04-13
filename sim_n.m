% [y,t,u] = sim_n(config, Ts)
function y = sim_n(config, Ts, nsim)
    
    config_ = config;
    y   = [];
    for i = 1:nsim
        % calculo `ek`
        % ek  = x0 - x_target;

        % calculo comando `dtk`
        % dtk = mpc_dualmode_switching(ek,H,Hf,Phi1Np,Qbar,Rbar,Lbar,cbar,Pf,Sf,bf,PhiNp,p);

        % simulando dinamica
        [y_,~,~]  = sim_1(config_, Ts);
        config_.x0 = y_(end,:)';

        y   = [y;y_];
    end
end