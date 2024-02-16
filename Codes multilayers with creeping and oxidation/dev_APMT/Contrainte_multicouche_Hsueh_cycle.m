function [dist]=Contrainte_multicouche_Hsueh_cycle(x)

global t_s;%thickness
global t_preoxide;
global t_tgo_top;
global t_tgo_bottom;
global t_c_top;
global t_c_bottom;
global Y_substrate;
global Y_coating;
global Y_tgo;
global Y_oxide;
global growthstrain_lat_top;
global growthstrain_lat_bottom;
global Lateral_preoxide;
global eps_thermique_tgo_top;
global eps_thermique_tgo_bottom;
global eps_thermique_coating_top;%misfit CTE
global eps_thermique_coating_bottom;
global eps_thermique_preoxide_top;
global eps_thermique_preoxide_bottom;
global eps_thermique_substrate;
global int_creep_tgo_top;
global int_creep_tgo_bottom
global int_creep_coating_top;%creep strain
global int_creep_coating_bottom;
global int_creep_preoxide_top;
global int_creep_preoxide_bottom;
global int_creep_substrate;
global moment_eps_tgo_top;
global moment_eps_tgo_bottom;
global moment_eps_preoxide_top;
global moment_eps_preoxide_bottom;
global moment_eps_substrate;
global moment_eps_coating_top;
global moment_eps_coating_bottom;
global moment_creep_tgo_top;
global moment_creep_tgo_bottom;
global moment_creep_preoxide_top;
global moment_creep_preoxide_bottom;
global moment_creep_substrate;
global moment_creep_coating_top;
global moment_creep_coating_bottom;
global eps_tgo_top1;
global eps_tgo_bottom1;
global C;
global strain_AlN;
global j;

c=x(1);
t_b=x(2);
r=x(3);

%formalisme de Hsueh
F7=Y_tgo*(c*t_tgo_top(j)-eps_tgo_top1*t_tgo_top(j)-eps_thermique_tgo_top-int_creep_tgo_top-t_tgo_top(j)*growthstrain_lat_top(j));
F6=Y_coating*(c*t_c_top(j)-eps_thermique_coating_top-int_creep_coating_top-t_c_top(j)*strain_AlN-t_c_top(j)*C);
F5=Y_oxide*(c*t_preoxide-eps_thermique_preoxide_top-int_creep_preoxide_top-t_preoxide*Lateral_preoxide);
F4=Y_substrate*(t_s*c-eps_thermique_substrate-int_creep_substrate);
F3=Y_oxide*(c*t_preoxide-eps_thermique_preoxide_bottom-int_creep_preoxide_bottom-t_preoxide*Lateral_preoxide);
F2=Y_coating*(c*t_c_bottom(j)-eps_thermique_coating_bottom-int_creep_coating_bottom-t_c_bottom(j)*strain_AlN-t_c_bottom(j)*C);
F1=Y_tgo*(c*t_tgo_bottom(j)-eps_tgo_bottom1*t_tgo_bottom(j)-eps_thermique_tgo_bottom-int_creep_tgo_bottom-t_tgo_bottom(j)*growthstrain_lat_bottom(j));
sum_biax_forces=F1+F2+F3+F4+F5+F6+F7;

%formalisem de Hsueh
C1 = @(z)Y_tgo*(z-t_b)/r;
C2 = @(z)Y_coating*(z-t_b)/r;
C3 = @(z)Y_oxide*(z-t_b)/r;
C4 = @(z)Y_substrate*(z-t_b)/r;
C5 = @(z)Y_oxide*(z-t_b)/r;
C6 = @(z)Y_coating*(z-t_b)/r;
C7 = @(z)Y_tgo*(z-t_b)/r;
%sum_curv_forces=integral(C1,-(t_s+t_tgo(j)),-t_tgo(j))+integral(C2,-t_tgo(j),0)+integral(C3,0,t_c1)+integral(C4,t_c1,(t_c1+t_c2));
sum_curv_forces=integral(C1,-t_s/2-t_c_bottom(j)-t_tgo_bottom(j)-t_preoxide,-t_s/2-t_preoxide-t_c_bottom(j))+...
                integral(C2,-t_s/2-t_preoxide-t_c_bottom(j),-t_s/2)-t_preoxide+...
                integral(C3,-t_s/2-t_preoxide,-t_s/2)+...
                integral(C4,-t_s/2,t_s/2)+...
                integral(C5,t_s/2,t_s/2+t_preoxide)+...
                integral(C6,t_s/2+t_preoxide,t_s/2+t_preoxide+t_c_top(j))+...
                integral(C7,t_s/2+t_preoxide+t_c_top(j),t_s/2+t_preoxide+t_c_top(j)+t_tgo_top(j));

%formalisme de Hsueh
int_c_tgo_top=c/2*((t_tgo_top(j)+t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_c_top(j)+t_preoxide+t_s/2-t_b)^2);
int_c_coating_top=c/2*((t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_s/2+t_preoxide-t_b)^2);
int_c_preoxide_top=c/2*((t_s/2+t_preoxide-t_b)^2-(t_s/2-t_b)^2);
int_c_substrate=c/2*((t_s/2-t_b)^2-(-t_s/2-t_b)^2);
int_c_preoxide_bottom=c/2*((-t_s/2-t_b)^2-(-t_s/2-t_preoxide-t_b)^2);
int_c_coating_bottom=c/2*((-t_s/2-t_preoxide-t_b)^2-(-t_s/2-t_preoxide-t_c_bottom(j)-t_b)^2);
int_c_tgo_bottom=c/2*((-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2-(-t_tgo_bottom(j)-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2);

int_b_tgo_top=1/(3*r)*((t_tgo_top(j)+t_c_top(j)+t_preoxide+t_s/2-t_b)^3-(t_c_top(j)+t_preoxide+t_s/2-t_b)^3);
int_b_coating_top=1/(3*r)*((t_c_top(j)+t_preoxide+t_s/2-t_b)^3-(t_s/2+t_preoxide-t_b)^3);
int_b_preoxide_top=1/(3*r)*((t_s/2+t_preoxide-t_b)^3-(t_s/2-t_b)^3);
int_b_substrate=1/(3*r)*((t_s/2-t_b)^3-(-t_s/2-t_b)^3);
int_b_preoxide_bottom=1/(3*r)*((-t_s/2-t_b)^3-(-t_s/2-t_preoxide-t_b)^3);
int_b_coating_bottom=1/(3*r)*((-t_s/2-t_preoxide-t_b)^3-(-t_s/2-t_preoxide-t_c_bottom(j)-t_b)^3);
int_b_tgo_bottom=1/(3*r)*((-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^3-(-t_tgo_bottom(j)-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^3);

int_elaboration_AlN_top=strain_AlN/2*((t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_s/2+t_preoxide-t_b)^2);
int_elaboration_AlN_bottom=strain_AlN/2*((-t_s/2-t_preoxide-t_b)^2-(-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2);

int_preoxidation_AlN_top=C/2*((t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_s/2+t_preoxide-t_b)^2);
int_preoxidation_AlN_bottom=C/2*((-t_s/2-t_preoxide-t_b)^2-(-t_s/2-t_preoxide-t_c_bottom(j)-t_b)^2);

int_g_tgo_top=growthstrain_lat_top(j)/2*((t_tgo_top(j)+t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_c_top(j)+t_preoxide+t_s/2-t_b)^2);
int_g_tgo_bottom=growthstrain_lat_bottom(j)/2*((-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2-(-t_tgo_bottom(j)-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2);

int_eps_tgo_top=eps_tgo_top1/2*((t_tgo_top(j)+t_c_top(j)+t_preoxide+t_s/2-t_b)^2-(t_c_top(j)+t_preoxide+t_s/2-t_b)^2);
int_eps_tgo_bottom=eps_tgo_bottom1/2*((-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2-(-t_tgo_bottom(j)-t_c_bottom(j)-t_preoxide-t_s/2-t_b)^2);

int_g_preoxidation_top=Lateral_preoxide/2*((t_s/2+t_preoxide-t_b)^2-(t_s/2-t_b)^2);
int_g_preoxidation_bottom=Lateral_preoxide/2*((-t_s/2-t_b)^2-(-t_s/2-t_preoxide-t_b)^2);

M1=Y_tgo*...
   (int_c_tgo_bottom+int_b_tgo_bottom-...
   (moment_eps_tgo_bottom-t_b*eps_thermique_tgo_bottom)-...
   (moment_creep_tgo_bottom-t_b*int_creep_tgo_bottom)-...
   int_g_tgo_bottom-...
   int_eps_tgo_bottom);
M2=Y_coating*...
   (int_c_coating_bottom+int_b_coating_bottom-...
   (moment_eps_coating_bottom-t_b*eps_thermique_coating_bottom)-...
   (moment_creep_coating_bottom-t_b*int_creep_coating_bottom)-...
   int_elaboration_AlN_bottom-...
   int_preoxidation_AlN_bottom);
M3=Y_oxide*...
   (int_c_preoxide_bottom+int_b_preoxide_bottom-...
   (moment_eps_preoxide_bottom-t_b*eps_thermique_preoxide_bottom)-...
   (moment_creep_preoxide_bottom-t_b*int_creep_preoxide_bottom)-...
   int_g_preoxidation_bottom);
M4=Y_substrate*...
   (int_c_substrate+int_b_substrate-...
   (moment_eps_substrate-t_b*eps_thermique_substrate)-...
   (moment_creep_substrate-t_b*int_creep_substrate));
M5=Y_oxide*...
   (int_c_preoxide_top+int_b_preoxide_top-...
   (moment_eps_preoxide_top-t_b*eps_thermique_preoxide_top)-...
   (moment_creep_preoxide_top-t_b*int_creep_preoxide_top)-...
   int_g_preoxidation_top);
M6=Y_coating*...
   (int_c_coating_top+int_b_coating_top-...
   (moment_eps_coating_top-t_b*eps_thermique_coating_top)-...
   (moment_creep_coating_top-t_b*int_creep_coating_top)-...
   int_elaboration_AlN_top-...
   int_preoxidation_AlN_top);
M7=Y_tgo*...
   (int_c_tgo_top+int_b_tgo_top-...
   (moment_eps_tgo_top-t_b*eps_thermique_tgo_top)-...
   (moment_creep_tgo_top-t_b*int_creep_tgo_top)-...
   int_g_tgo_top-...
   int_eps_tgo_top);

sum_curv_moments=M1+M2+M3+M4+M5+M6+M7;
dist=abs(sum_biax_forces)+1000*abs(sum_curv_forces)+1e6*abs(sum_curv_moments);