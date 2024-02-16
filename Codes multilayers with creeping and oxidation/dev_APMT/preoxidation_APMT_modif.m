%Calculation of stress in a multilayer coating
function [C_preoxidation,c,t_b,r,t_s,t_preoxide,Creep_preoxidation_substrate,Creep_preoxidation_oxide,Lateral_preoxide]=preoxidation_APMT_modif(t_s_int)

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
global alpha_coating;
global alpha_substrate;
global alpha_oxide;
global Y_oxide;
global Y_substrate;
global Y_coating;
global T_ambiant;
global T_elaboration;
global T_preoxidation;
global lateral_strain_constant_preoxidation;
holdtime=24;%h
delta_t=1/120;%h
nb_t=holdtime/delta_t+1;
NA=6.022e23;
R=8.314472;
T=T_preoxidation;%°C, oxidation temperature

t_sub=t_s_int;

%physical proprities of oxide
M_oxide=101.96;%g/mol
density_oxide=3.95e6;%g/m3
atomic_volume_oxide=M_oxide/(2*NA*density_oxide);

%creep proprities of oxide
%B_oxide=exp(15.68301);%MPa-n.s-1
%n_oxide=1.64812;
%Q_oxide=460e3;
B_oxide=3e-9;%MPa-n.s-1
n_oxide=1;
Q_oxide=0;

%physical proptities of substrate
M_metal=26.98;%g/mole
density_substrate=7.25e6;%g/m3
mass_fraction_metal=0.05;
atomic_volume_metal=M_metal/(NA*density_substrate*mass_fraction_metal);

%creep proprities of substrate
%B_substrate=5.96e6;%MPa-n.s-1
%n_substrate=5.5;
%Q_substrate=392e3;%J/mol
B_substrate=78.978;%MPa-n.s-1, MA956
n_substrate=4.9827;
Q_substrate=453e3;%J/mol

%oxidation proprities
kp=2.366e-17;%m2.s-1

stress_s=0;
creep_s=0;
creep_rate_s=0;
th_s(1)=t_sub;
creep_o=0;
creep_rate_o=0;
lateral_strain=0;
time=0;
radis=inf;
tb=0;

for j=2:nb_t
    time(j)=time(j-1)+delta_t;
    %update oxide scale thickness
    t_oxide(j)=(kp*time(j)*3600)^(1/2);%m
    %update substrate thickness
    d_oxide=t_oxide(j)-t_oxide(j-1);%m
    d_substrate=-d_oxide*(atomic_volume_metal/atomic_volume_oxide)^(1/3);%m
    th_s(j)=th_s(j-1)+d_substrate;

    %update creep strain for oxide and substrate
    creep_o(j)=creep_rate_o(j-1)*delta_t+creep_o(j-1);
    creep_s(j)=creep_rate_s(j-1)*delta_t+creep_s(j-1);

    %update lateral growth strain for oxide
    lateral_strain_rate=1/2*lateral_strain_constant_preoxidation*kp^(1/2)*(3600*time(j))^(-1/2);
    lateral_strain(j)=lateral_strain(j-1)+lateral_strain_rate*delta_t*3600;

    %sum of biaxial stress=0
    sum_biaxial_stress= @(cst) 2*Y_oxide*(cst-creep_o(j)-lateral_strain(j))*t_oxide(j)+...
        Y_substrate*(cst-creep_s(j))*th_s(j);
    cst=fzero(sum_biaxial_stress,0);

    %stress calculation
    stress_s(j)=Y_substrate*(cst-creep_s(j));
    stress_o(j)=Y_oxide*(cst-creep_o(j)-lateral_strain(j));
    %update creep rate of oxide and substrate
    creep_rate_s(j)=3600*sign(stress_s(j))*B_substrate*(abs(stress_s(j)/1e6))^n_substrate*exp(-Q_substrate/(R*(273.15+T)));
    creep_rate_o(j)=3600*sign(stress_o(j))*B_oxide*(abs(stress_o(j))/1e6)^n_oxide*exp(-Q_oxide/(R*(273.15+T)));

end

% figure(1)
% plot(time,stress_o);
% figure(2)
% plot(time,stress_s);
% figure(3)
% plot(time,creep_rate_s);
% figure(4)
% plot(time,creep_rate_o);
% figure(5)
% plot(time,lateral_strain);
% figure(6)
% plot(time,creep_o);
% figure(7)
% plot(time,creep_s);

t_sub=th_s(j);
t_oxide=t_oxide(j);
Creep_oxide=creep_o(j);
Creep_substrate=creep_s(j);
Lateral_preoxidation=lateral_strain(j);
eps_preoxide=integral(alpha_oxide,T,T_elaboration);
eps_substrate=integral(alpha_substrate,T,T_elaboration);

%sum of biaxial stress=0
sum_biaxial_stress= @(cst) 2*Y_oxide*(cst-Creep_oxide-Lateral_preoxidation-eps_preoxide)*t_oxide+...
    Y_substrate*(cst-Creep_substrate-eps_substrate)*t_sub;%+...
%Y_coating*(cst-strain_AlN)*t_c1+...
%Y_coating*(cst-strain_AlN)*t_c2;
cst=fzero(sum_biaxial_stress,0);

disp('Initial configuration')
disp(['c=',num2str(cst),' tb=',num2str(tb),' R=',num2str(radis)]);

C=cst;

if (t_c1-t_c2)==0
    %sum of biaxial stress=0
    sum_biaxial_stress= @(cst) 2*Y_oxide*(cst-Creep_oxide-Lateral_preoxidation-eps_preoxide)*t_oxide+...
        Y_substrate*(cst-Creep_substrate-eps_substrate)*t_sub+...
        2*Y_coating*(cst-strain_AlN-C)*t_c1;
    cst=fzero(sum_biaxial_stress,0);
    disp(['c=',num2str(cst),' tb=',num2str(tb),' R=',num2str(radis)]);
end

eps_coating=integral(alpha_coating,T_elaboration,T_ambiant);
eps_preoxide=integral(alpha_oxide,T_preoxidation,T_ambiant);
eps_substrate=integral(alpha_substrate,T_preoxidation,T_ambiant);

if (t_c1-t_c2)==0
    %sum of biaxial stress=0
    sum_biaxial_stress= @(c) 2*Y_oxide*(c-Creep_oxide-Lateral_preoxidation-eps_preoxide)*t_oxide+...
        Y_substrate*(c-Creep_substrate-eps_substrate)*t_sub+...
        2*Y_coating*(c-eps_coating-strain_AlN-C)*t_c1;
    c=fzero(sum_biaxial_stress,0);
    t_b=0;
    r=Inf;
    disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);

else

    fval=0;
    while fval==0
        [x,fval]=GA(@Contrainte_multicouche_Hsueh_preoxidation,3,1e-7);
    end

    disp('Converging, next step beginning');
    c=x(1);
    t_b=x(2);
    r=x(3);
    % disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);
    % if x(3)>0
    %     disp('Dépot en compression');
    % end %definition de Hsueh: compression r>0
    % if x(3)<0
    %     disp('Dépot en tension');
    % end
    % disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);
end
t_s=t_sub;
t_preoxide=t_oxide;
C_preoxidation=C;
Creep_preoxidation_oxide=Creep_oxide;
Creep_preoxidation_substrate=Creep_substrate;
Lateral_preoxide=Lateral_preoxidation;

%stress_s_amb=Y_substrate*(c+(-t_b)/r-creep_s(j)-eps_substrate);
%stress_o_amb=Y_oxide*(c+(t_oxide+t_sub/2-t_b)/r-creep_o(j)-lateral_strain(j)-eps_preoxide);
%stress_c_amb=Y_coating*(c+(t_c1+t_oxide+t_sub/2-t_b)/r-eps_coating-strain_AlN-C);
