clc
clear
close all
rng('shuffle', 'twister')
warning('off','all')
format 'ShortE'
global t_s;%thickness
global t_c1;
global t_c2;
global t_preoxide;
global strain_AlN;
global C;
global Creep_preoxidation_oxide;
global Creep_preoxidation_substrate;
global Lateral_preoxide;
global alpha_coating;
global alpha_substrate;
global alpha_tgo;
global alpha_oxide;
global Y_oxide;
global Y_substrate;
global Y_coating;
global Y_tgo;
global t_int;
global nb_c1;
global nb_c2;
global nb_s;
global nb_o;
global nb_preoxide;
global nb;
global nb_t;
global deltat;
global T_preoxidation;
global T_elaboration;
global T_bottom;
global T_ambiant;
global q_thermique;
global size_z_sub;
global size_z_preoxide;
global size_z_coating;
global size_z_tgo;
global lateral_strain_constant_preoxidation;
global lateral_strain_constant_AlN;
global k;

t_s_int=3e-3;%m
t_c1=21e-6;
t_c2=27e-6;
strain_AlN=-3e-3;
T_preoxidation=1100;
T_elaboration=1100;%°C température pour la croissance d'AlN
T_ambiant=25;%°C température finale d'une cycle
%T_top=1400-273.15;%°C température fixée pour le top du système
T_bottom=1000;%°C
%T_int=[T_top T_ambiant];
q_thermique=-1000e3;%W.m-2

t_int(1)=0;%tgo_top
t_int(2)=0;%tgo_bottom
t_int(3)=t_c1;%c_top
t_int(4)=t_c2;%c_bottom

holdtime=1;%h
deltat=1/10;%h
nb_t=3+holdtime/deltat;
nb_cycle=10;

%Data for substrate, APMT
E_substrate=220e9;%Pa, ambient temperature
v_substrate=0.3;
Y_substrate=E_substrate/(1-v_substrate);
a0_s=1.07677e+01;%Kanthal APMT datasheet and Handbook APMT, polyfit by order 5
a1_s=4.05923e-03;
a2_s=-5.85985e-08;
a3_s=-6.17356e-09;
a4_s=1.06511e-11;
a5_s=-4.56942e-15;

%Data for coating, AlN
E_coating=340e9;%Pa
nu_coating=0.21;%http://accuratus.com/silicar.html
Y_coating=E_coating/(1-nu_coating);
a0_c=2.47411e+00;
a1_c=1.11737e-02;
a2_c=-1.32814e-05;
a3_c=8.67958e-09;
a4_c=-2.94971e-12;
a5_c=4.05708e-16;

%Data for oxide, Al2O3
E_oxide=379e9;%Pa, X.Dong 2012 oxidation stress evolution and relation of oxide film/metal substrate system
v_oxide=0.25;
Y_oxide=E_oxide/(1-v_oxide);
a0_oxide=5.07720e+00;
a1_oxide=1.43166e-02;
a2_oxide=-2.30537e-05;
a3_oxide=2.06357e-08;
a4_oxide=-8.83916e-12;
a5_oxide=1.42469e-15;

Y_tgo=Y_oxide;
alpha_substrate = @(T)(a0_s+a1_s.*T+a2_s.*T.^2+a3_s.*T.^3+a4_s.*T.^4+a5_s.*T.^5)/1000000;%polyfit by order 5
alpha_oxide = @(T)(a0_oxide+a1_oxide.*T+a2_oxide.*T.^2+a3_oxide.*T.^3+a4_oxide.*T.^4+a5_oxide.*T.^5)/1000000;%polyfit by order 5
alpha_coating = @(T)(a0_c+a1_c.*T+a2_c.*T.^2+a3_c.*T.^3+a4_c.*T.^4+a5_c.*T.^5)/1000000;%polyfit by order 5
alpha_tgo = alpha_oxide;

eps_coating_depot=integral(alpha_coating,T_elaboration,T_ambiant);
eps_preoxide_depot=integral(alpha_oxide,T_preoxidation,T_ambiant);
eps_substrate_depot=integral(alpha_substrate,T_preoxidation,T_ambiant);

lateral_strain_constant_preoxidation=4.25e4;%m-1
lateral_strain_constant_AlN=0.5e-4;%m-1

[C_preoxidation,c,t_b,r,t_s,t_preoxide,Creep_preoxidation_substrate,Creep_preoxidation_oxide,Lateral_preoxide]=preoxidation_APMT_modif(t_s_int);

% Y_coating*(c+(t_c1+t_preoxide+t_s/2-t_b)/r-eps_coating_depot-strain_AlN-C);
% Y_substrate*(c-t_b/r-eps_substrate_depot-Creep_preoxidation_substrate);

size_z_sub=2e-8;
nb_s=ceil(t_s/size_z_sub);
size_z_coating=1e-7;
nb_c1=ceil(t_c1/size_z_coating)-1;
nb_c2=ceil(t_c2/size_z_coating)-1;
size_z_preoxide=1e-8;
nb_preoxide=ceil(t_preoxide/size_z_preoxide);
size_z_tgo=1e-8;
nb_o=1000;
nb=1+2*nb_o+nb_c1+nb_c2+nb_s+2*nb_preoxide;

creepstrain_int=zeros(nb,1);

for i=1:nb%z position
    z(1:nb_o+1,1)=t_s/2+t_c1+t_preoxide;
    if (nb_o+1)<i&&i<=(nb_c1+nb_o+1)%c1
        z(i,1)=t_s/2+t_c1+t_preoxide-size_z_coating*(i-nb_o-1);
        stress_elaboration(i,1)=Y_coating*(c+(z(i,1)-t_b)/r-eps_coating_depot-strain_AlN-C);
    end
    if (nb_o+nb_c1+1)<i&&i<=(nb_c1+nb_o+nb_preoxide+1)%preoxide_top
        z(i,1)=t_s/2-size_z_preoxide*(i-nb_o-nb_c1-nb_preoxide-1);
        stress_elaboration(i,1)=Y_oxide*(c+(z(i,1)-t_b)/r-eps_preoxide_depot-Creep_preoxidation_oxide-Lateral_preoxide);
        creepstrain_int(i,1)=Creep_preoxidation_oxide;
    end
    if (nb_preoxide+nb_c1+nb_o+1)<i&&i<(nb_s+nb_preoxide+nb_c1+nb_o+1)%substrate
        z(i,1)=-size_z_sub*(i-nb_s/2-1-nb_c1-nb_o-nb_preoxide);
        stress_elaboration(i,1)=Y_substrate*(c+(z(i,1)-t_b)/r-eps_substrate_depot-Creep_preoxidation_substrate);
        creepstrain_int(i,1)=Creep_preoxidation_substrate;
    end
    %z(nb_s+nb_preoxide+nb_c1+nb_o+1,1)=-t_s/2;
    %stress_elaboration(nb_s+nb_preoxide+nb_c1+nb_o+1,1)=Y_oxide*(c+(z(nb_s+nb_preoxide+nb_c1+nb_o+1,1)-t_b)/r-eps_preoxide_depot-Creep_preoxidation_oxide-Lateral_preoxide);
    %creepstrain_int(nb_s+nb_preoxide+nb_c1+nb_o+1,1)=Creep_preoxidation_oxide;
    if (nb_s+nb_preoxide+nb_c1+nb_o+1)<=i&&i<(nb_s+2*nb_preoxide+nb_c1+nb_o+1)%preoxide_bottom
        z(i,1)=-t_s/2-size_z_preoxide*(i-nb_s-nb_c1-nb_o-nb_preoxide-1);
        stress_elaboration(i,1)=Y_oxide*(c+(z(i,1)-t_b)/r-eps_preoxide_depot-Creep_preoxidation_oxide-Lateral_preoxide);
        creepstrain_int(i,1)=Creep_preoxidation_oxide;
    end
    %z(nb_s+2*nb_preoxide+nb_c1+nb_o+1,1)=-t_s/2-t_preoxide;
    if (nb_s+2*nb_preoxide+nb_c1+nb_o+1)<=i&&i<(nb_s+2*nb_preoxide+nb_c1+nb_c2+nb_o+1)%c2
        z(i,1)=-t_s/2-t_preoxide-size_z_coating*(i-nb_s-2*nb_preoxide-nb_c1-nb_o-1);
        stress_elaboration(i,1)=Y_coating*(c+(z(i,1)-t_b)/r-eps_coating_depot-strain_AlN-C);
    end
    z((nb-nb_o):nb,1)=-t_s/2-t_c2-t_preoxide;
    stress_elaboration((nb-nb_o+1):nb,1)=0;
end

for k=1:nb_cycle
    if k==1
        T_oxi_int=0;
        z_int=z;
        tao_int=0;
        growthstrain_lat_top_int=0;
        growthstrain_lat_bottom_int=0;
        th_int=t_int';
        eps_tgo_top_int=0;
        eps_tgo_bottom_int=0;
        stress_int=stress_elaboration;
        stress_ambiant=stress_elaboration;
        time_ambiant=0;
        %creepstrain_int=zeros(nb,1);
    end

    [T,T_oxi_out,z_out,tao_accu,th_cycle,growthstrain_lat_top_accu,growthstrain_lat_bottom_accu,eps_tgo_top_out,eps_tgo_bottom_out,stress_cycle,creepstrain_cycle,...
        Curvature_cycle]=cycle_GA(T_oxi_int,z_int,tao_int,th_int,growthstrain_lat_top_int,growthstrain_lat_bottom_int,eps_tgo_top_int,eps_tgo_bottom_int,stress_int,...
        creepstrain_int);

    T_oxi_int=T_oxi_out;
    z_int=z_out;
    tao_int=tao_accu;
    growthstrain_lat_top_int=growthstrain_lat_top_accu;
    growthstrain_lat_bottom_int=growthstrain_lat_bottom_accu;
    th_int(1,1)=th_cycle(1,nb_t);
    th_int(2,1)=th_cycle(2,nb_t);
    th_int(3,1)=th_cycle(3,nb_t);
    th_int(4,1)=th_cycle(4,nb_t);
    stress_int=stress_cycle(:,nb_t);
    creepstrain_int=creepstrain_cycle(:,nb_t);
    eps_tgo_top_int=eps_tgo_top_out;
    eps_tgo_bottom_int=eps_tgo_bottom_out;

    for j=2:(nb_t-1)
        temps(j)=1/12+deltat*(j-2);
    end
    temps(nb_t)=temps(nb_t-1)+1/12;

    if k==1
        Temperature=T;
        stress=stress_cycle;
        creepstrain=creepstrain_cycle;
        Curvature_cycle(1)=1/r;
        Curvature=Curvature_cycle;
        thickness=th_cycle;
        TIME=temps;
        stress_ambiant=[stress_ambiant,stress_cycle(:,nb_t)];
        time_ambiant=[time_ambiant,time_ambiant(k)+TIME(nb_t)];
    end

    if k>1
        thickness=[thickness,th_cycle];
        Temperature=[Temperature,T];
        stress=[stress,stress_cycle];
        creepstrain=[creepstrain,creepstrain_cycle];
        Curvature_cycle(1)=Curvature(nb_t*(k-1));
        Curvature=[Curvature,Curvature_cycle];
        TIME=[TIME,temps+TIME(nb_t*(k-1))];
        stress_ambiant=[stress_ambiant,stress_cycle(:,nb_t)];
        time_ambiant=[time_ambiant,time_ambiant(k)+TIME(nb_t)];
    end
    disp(['****************************Cycle ',num2str(k), '***************' ...
        '' ...
        '*************'])
end

plot(time_ambiant,stress_ambiant(1002,:))
drawnow
fid=fopen('APMT_multilayer coating_oxidation_Curvature_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
fprintf(fid,'%12.5f\n',Curvature);
fclose(fid);

fid=fopen('APMT_multilayer coating_oxidation_TIME_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
fprintf(fid,'%12.5f\n',TIME);
fclose(fid);

fid=fopen('APMT_multilayer coating_oxidation_Timeambiant_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
fprintf(fid,'%12.5f\n',time_ambiant);
fclose(fid);

[m,n]=size(stress);
fid=fopen('APMT_multilayer coating_oxidation_stress_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
for i=1:1:m
    for j=1:1:n
        if j==n
            fprintf(fid,'%12.5f\n',stress(i,j));
        else
            fprintf(fid,'%12.5f\t',stress(i,j));
        end
    end
end
fclose(fid);

[m,n]=size(creepstrain);
fid=fopen('APMT_multilayer coating_oxidation_creepstrain_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
for i=1:1:m
    for j=1:1:n
        if j==n
            fprintf(fid,'%12.5e\n',creepstrain(i,j));
        else
            fprintf(fid,'%12.5e\t',creepstrain(i,j));
        end
    end
end
fclose(fid);

[m,n]=size(stress_ambiant);
fid=fopen('APMT_multilayer coating_oxidation_stressambiant_Tuniforme_maille50000_corr_linear oxidation_n(AlN)1_creep(AlN)times200.txt','wt');
for i=1:1:m
    for j=1:1:n
        if j==n
            fprintf(fid,'%12.5f\n',stress_ambiant(i,j));
        else
            fprintf(fid,'%12.5f\t',stress_ambiant(i,j));
        end
    end
end
fclose(fid);
