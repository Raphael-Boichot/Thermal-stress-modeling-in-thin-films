function [dist]=Contrainte_multicouche_Hsueh_preoxidation(x)
global t_sub;
global t_c1;
global t_c2;
global t_oxide;
global strain_AlN;
global C;
global Creep_oxide;
global Creep_substrate;
global Lateral_preoxidation;
global eps_coating;
global eps_preoxide;
global eps_substrate;
global Y_oxide;
global Y_substrate;
global Y_coating;

c=x(1);
t_b=x(2);
r=x(3);

%formalisme de Hsueh
F1=Y_substrate*(c-eps_substrate-Creep_substrate)*t_sub;
F2=Y_coating*(c-eps_coating-strain_AlN-C)*t_c1;
F3=Y_coating*(c-eps_coating-strain_AlN-C)*t_c2;
F4=Y_oxide*(c-eps_preoxide-Creep_oxide-Lateral_preoxidation)*t_oxide;
sum_biax_forces=F1+F2+F3+2*F4;

C1 = @(z)Y_substrate*(z-t_b)/r;
C2 = @(z)Y_coating*(z-t_b)/r;
C3 = @(z)Y_oxide*(z-t_b)/r;
sum_curv_forces=integral(C1,-t_sub/2,t_sub/2)+integral(C2,(t_sub/2+t_oxide),(t_c1+t_sub/2+t_oxide))+...
                integral(C2,-(t_sub/2+t_c2+t_oxide),(-t_sub/2+t_oxide))+integral(C3,-(t_sub/2+t_oxide),-t_sub/2)+...
                integral(C3,t_sub/2,(t_sub/2+t_oxide));

%formalisme de Hsueh
M1 = @(z)Y_substrate*(c+(z-t_b)/r-eps_substrate-Creep_substrate).*(z-t_b);
M2 = @(z)Y_coating*(c+(z-t_b)/r-eps_coating-strain_AlN-C).*(z-t_b);
M3 = @(z)Y_oxide*(c+(z-t_b)/r-eps_preoxide-Creep_oxide-Lateral_preoxidation).*(z-t_b);
sum_curv_moments=integral(M1,-t_sub/2,t_sub/2)+integral(M2,(t_sub/2+t_oxide),(t_c1+t_sub/2+t_oxide))+...
                 integral(M2,-(t_sub/2+t_c2+t_oxide),(-t_sub/2+t_oxide))+integral(M3,-(t_sub/2+t_oxide),-t_sub/2)+...
                 integral(M3,t_sub/2,(t_sub/2+t_oxide));

dist=abs(sum_biax_forces)+1000*abs(sum_curv_forces)+1e6*abs(sum_curv_moments);