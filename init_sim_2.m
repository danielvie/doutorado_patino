%%
% --------------------------------------
% CONFIG
% --------------------------------------
%{
	Mode	u1 	u2 	u3
	0		0  	0 	0
	1		0  	0 	1
	2		0  	1 	0
	3		0  	1 	1
	4		1  	0 	0
	5		1  	0 	1
	6		1  	1 	0
	7		1  	1 	1
%}

function config = init_sim_2()
	% equacoes de estado Buck-Boost converter
	R  = 10; % Ω
	L  = 10; % mH
	C1 = 40; % μF
	C2 = 40; % μF
	E  = 30; % V

	A0  = [
		0, 0, 0
        0, 0, 0
        0, 0, -R/L
	];
	
	A1  = [
		0, 0, 0
        0, 0, 1/C2
        0, -1/L, -R/L
	];
	
	A2  = [
		0, 0, 1/C1
        0, 0, 1/C2
        -1/L, 1/L, -R/L
	];

	A3  = [
		0, 0, 1/C1
        0, 0, 0
        -1/L, 1/L, -R/L
	];

	A4  = [
		0, 0, -1/C1
        0, 0, 0
        1/L, 0, -R/L
	];

	A5  = [
		0, 0, -1/C1
        0, 0, 1/C2
        1/L, -1/L, -R/L
	];

	A6  = [
		0, 0, 0
        0, 0, -1/C2
        0, 1/L, -R/L
	];

	A7  = [
		0, 0, 0
        0, 0, 0
        0, 0, -R/L
	];

	B0 = [0; 0; E/L];

	B1 = [0; 0; 0];

	B2 = [0; 0; E/L];

	B3 = [0; 0; 0];

	B4 = [0; 0; E/L];

	B5 = [0; 0; 0];

	B6 = [0; 0; E/L];

	B7 = [0; 0; 0];
		
	Cc = eye(3);

	Dc = [0;0;0];

	config.modo  = [0, 1, 2, 3, 4, 5, 6, 7];
	config.Ac    = {A0, A1, A2, A3, A4, A5, A6, A7};
	config.Bc    = {B0, B1, B2, B3, B4, B5, B6, B7};
	config.Cc    = Cc;
	config.Dc    = Dc;

	config.tstep = 0.0001;
	config.xref  = [2/3*E, 1/3*E, 1];
	
	config.Q     = diag([10, 5, 20000]);
	config.Tpmax = 0.0004; % 0.4ms
    config.tmin  = 0.000022; % 0.022 ms
    
	config.smax  = 12;
	config.x0    = [9.9247; 19.2928; 0.9823];
			
end