function [dist]=Contrainte_multicouche_Hsueh_4couches(x);


global t_c;%thickness
global t_s;
global eps_thermique_coating_top;%misfit CTE
global eps_thermique_coating_bottom;
global eps_thermique_substrate;
global Y_substrate;
global Y_coating;
global int_creep_coating_top;%creep strain
global int_creep_coating_bottom;
global int_creep_substrate;
global moment_eps_substrate;
global moment_eps_coating_top;
global moment_eps_coating_bottom;
global moment_creep_substrate;
global moment_creep_coating_top;
global moment_creep_coating_bottom;
global j;
global strain_AlN;

c=x(1);
t_b=x(2);
r=x(3);

%formalisme de Hsueh

F2=Y_substrate*(t_s*c-eps_thermique_substrate-int_creep_substrate);
F3=Y_coating*(c*t_c-eps_thermique_coating_top-int_creep_coating_top-t_c*strain_AlN);
F1=Y_coating*(c*t_c-eps_thermique_coating_bottom-int_creep_coating_bottom-t_c*strain_AlN);
sum_biax_forces=F1+F2+F3;


%formalisem de Hsueh
C2 = @(z)Y_substrate*(z-t_b)/r;
C1 = @(z)Y_coating*(z-t_b)/r;
C3 = @(z)Y_coating*(z-t_b)/r;
%sum_curv_forces=integral(C1,-(t_s+t_tgo(j)),-t_tgo(j))+integral(C2,-t_tgo(j),0)+integral(C3,0,t_c1)+integral(C4,t_c1,(t_c1+t_c2));
sum_curv_forces=integral(C1,-t_s/2-t_c,-t_s/2)+integral(C2,-t_s/2,t_s/2)+integral(C3,t_s/2,t_s/2+t_c);


%formalisme de Hsueh
int_c_substrate=c/2*((t_s/2-t_b)^2-(-t_s/2-t_b)^2);
int_c_coating_top=c/2*((t_c+t_s/2-t_b)^2-(t_s/2-t_b)^2);
int_c_coating_bottom=c/2*((-t_s/2-t_b)^2-(-t_s/2-t_c-t_b)^2);
int_b_substrate=1/(3*r)*((t_s/2-t_b)^3-(-t_s/2-t_b)^3);
int_b_coating_top=1/(3*r)*((t_c+t_s/2-t_b)^3-(t_s/2-t_b)^3);
int_b_coating_bottom=1/(3*r)*((-t_s/2-t_b)^3-(-t_s/2-t_c-t_b)^3);
int_elaboration_AlN_top=strain_AlN/2*((t_c+t_s/2-t_b)^2-(t_s/2-t_b)^2);
int_elaboration_AlN_bottom=strain_AlN/2*((-t_s/2-t_b)^2-(-t_c-t_s/2-t_b)^2);

M2=Y_substrate*(int_c_substrate+int_b_substrate-(moment_eps_substrate-t_b*eps_thermique_substrate)-(moment_creep_substrate-t_b*int_creep_substrate));
M1=Y_coating*(int_c_coating_bottom+int_b_coating_bottom-(moment_eps_coating_bottom-t_b*eps_thermique_coating_bottom)-(moment_creep_coating_bottom-t_b*int_creep_coating_bottom)-int_elaboration_AlN_bottom);
M3=Y_coating*(int_c_coating_top+int_b_coating_top-(moment_eps_coating_top-t_b*eps_thermique_coating_top)-(moment_creep_coating_top-t_b*int_creep_coating_top)-int_elaboration_AlN_top);
%sum_curv_moments=M1+M2+M3+M4;
sum_curv_moments=M1+M2+M3;

dist=abs(sum_biax_forces)+100*abs(sum_curv_forces)+1e6*abs(sum_curv_moments);
%dist=abs(sum_biax_forces)+1000*abs(sum_curv_forces)+1e6*abs(sum_curv_moments);

