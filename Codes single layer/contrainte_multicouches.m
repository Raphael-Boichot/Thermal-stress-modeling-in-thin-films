function [dist]=contrainte_multicouches(x)

global ilots;
global t_f;
global t_s;
global eps_thermique_film;
global eps_thermique_substrate;
global Y_substrate;
global Y_film;
global eps_misfit;
global eps_coalescence;

c=x(1);
t_b=x(2);
r=x(3);

%attention, changement du signe pour alpha par rapport à Hsueh
F1=Y_substrate*(c+eps_thermique_substrate)*t_s;
F2=Y_film*(c+eps_thermique_film-eps_misfit-eps_coalescence*ilots)*t_f;
sum_biax_forces=F1+F2;

C1 = @(z)Y_substrate*(z-t_b)/r;
C2 = @(z)Y_film*(z-t_b)/r;
sum_curv_forces=integral(C1,-t_s,0)+integral(C2,0,t_f);

%attention, changement du signe pour alpha par rapport à Hsueh
M1 = @(z)Y_substrate*(c+(z-t_b)/r+eps_thermique_substrate).*(z-t_b);
M2 = @(z)Y_film*(c+(z-t_b)/r+eps_thermique_film-eps_misfit-eps_coalescence*ilots).*(z-t_b);
sum_curv_moments=integral(M1,-t_s,0)+integral(M2,0,t_f);

dist=abs(sum_biax_forces)+100*abs(sum_curv_forces)+1e6*abs(sum_curv_moments);

