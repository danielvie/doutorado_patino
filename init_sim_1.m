%%
% --------------------------------------
% CONFIG
% --------------------------------------

function config = init_sim_1()
      % equacoes de estado Buck-Boost converter
      R = 1;
      L = 1;
      C = 1;
      E = 1;

      A1 = [0,0
            0, -1/(R*C)];
      
      A2 = [   0,    1/L
            -1/C, -1/(R*C)];

      B1 = [E/L; 0];

      B2 = [0; 0];

      Cc = eye(2);

      Dc = [0;0];

      tmin  = 0.25;
      xref  = [2, -1];
      % xref: [vc_ref; il_ref]

      Q     = diag([1,1]);
      Tpmax = 1;
      smax  = 2;
      r     = 1;
      % r: 

      tstep = 0.0001;

      config.A1 = A1;
      config.A2 = A2;
      config.B1 = B1;
      config.B2 = B2;
      config.Cc = Cc;
      config.Dc = Dc;

      config.modes  = [0, 1]; % modo de operacao
      config.ur     = [1, 0]; % controle associado com modo
      
      config.Ac = {A1, A2};
      config.Bc = {B1, B2};
      
      config.tstep = tstep;
      config.xref  = xref;
      config.tmin  = tmin;
      config.Q     = Q;
      config.Tpmax = Tpmax;
      config.smax  = smax;
      config.r     = r;
            
end