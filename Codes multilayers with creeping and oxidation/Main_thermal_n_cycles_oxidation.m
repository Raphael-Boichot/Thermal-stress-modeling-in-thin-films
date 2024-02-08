% ------------------------
% 1-code 'Main_thermal_n_cycles_oxidation' & 'cycle_GA' & 'Contrainte_multicouche_Hsueh_4' couche are used for calculating the stress evolution during thermal cycling
% 
% 2-code 'minimize_GA_TZM_1' is used for calculating initial stress distribution in the system (just after deposition)
% 
% 3-files 'AlN_TZM_AlN_q10000_xx_CTEofT' are obtained using instantaneous CTE for claculations
%   also, we used a deltat of 1/12 (h) and a node number of 5000 for the substrate during these calculations
% 
% 4-files 'AlN_TZM_AlN_q10000' are obtained using constant CTE
%   also, we used a deltat of 1/12 (h) and a node number of 5000 for the substrate during these calculations
% 
% 5-files 'AlN_TZM_AlN_maille50000_q10000_Tambxx_xx': use 25 °C or 20 °C as ambiant temperature in the code (for the others without Tamb mentioned in the file name, we use 25 °C as ambiant temperature)
%   also, we used a deltat of 1/6 (h) and a node number of 50000 for the substrate during these calculations
% 
% 5-files 'AlN_TZM_AlN_maille50000_q10000_xx_CTEofT': use instantaneous CTE for claculations
%   also, we used a deltat of 1/6 (h) and a node number of 50000 for the substrate during these calculations
% ------------------------

clc
clear all
close all
rng('shuffle', 'twister')
format 'ShortE'

global t_c;%thickness
global t_s;
global nb_c;
global nb_s;
global nb;
global nb_t;
global deltat;
global strain_AlN;
global T_elaboration;
global q_thermique;

T_elaboration=1200;%°C température pour la croissance d'AlN
T_ambiant=20;%°C température finale d'une cycle
%T_top=1400-273.15;%°C température fixée pour le top du système
T_bottom=1000;%°C température fixée pour la bas du système
T_int=[T_bottom T_ambiant];
q_thermique=-10000e3;%W.m-2

%strain_AlN=-3e-3;%déformation de croissance
strain_AlN=0;

t_s=1e-3; %m, épaisseur du substrat
t_c=60e-6;

holdtime=10;%h
deltat=1/12;%h

nb_c=300;
nb_s=5000;
nb=1+2*nb_c+nb_s;

nb_t=3+holdtime/deltat;
%nb_t=4+holdtime/deltat;

nb_cycle=10;

[r,stress_elaboration]=minimize_GA_TZM_1(T_elaboration,T_ambiant);%initial stress just after deposition


for k=1:nb_cycle
    if k==1
        tao_int=0;
        %growthstrain_thick_int=0;
        %growthstrain_lat_int=0;
        %t_tgo_int=t_tgo;
        stress_int=stress_elaboration;
        creepstrain_int=zeros(nb,1);
    end
    
    [tao_accu,stress_cycle,creepstrain_cycle,Curvature_cycle]=cycle_GA(T_int,tao_int,stress_int,creepstrain_int);
    
    tao_int=tao_accu;
    stress_int=stress_cycle(:,nb_t);
    creepstrain_int=creepstrain_cycle(:,nb_t);
    
    for j=2:(nb_t-1)
        temps(j)=1/12+deltat*(j-2);
    end
    temps(nb_t)=temps(nb_t-1)+1/12;
    
    if k==1
        stress=stress_cycle;
        creepstrain=creepstrain_cycle;
        Curvature_cycle(1)=1/r;
        Curvature=Curvature_cycle;
        %th_tgo=th_tgo_cycle;
        TIME=temps;
    end
    
    if k>1
    
    %growthstrain_thick_int=growthstrain_thick_accu;
    %growthstrain_lat_int=growthstrain_lat_accu;
    %t_tgo_int=th_tgo_cycle(nb_t);
    
    %th_tgo=[th_tgo,th_tgo_cycle];
    stress=[stress,stress_cycle];
    creepstrain=[creepstrain,creepstrain_cycle];
    Curvature_cycle(1)=Curvature(nb_t*(k-1));
    Curvature=[Curvature,Curvature_cycle];  
    TIME=[TIME,temps+TIME(nb_t*(k-1))];
    end
    k
    
    figure(1)
    plot(TIME,Curvature,'LineWidth',1)
    figure(2)
    plot(TIME,stress(nb_c+1,:),'LineWidth',1)
    figure(3)
    plot(TIME,stress(nb_c+nb_s+1,:),'LineWidth',1)
    xlabel('TIME (h)')
    ylabel('Stress (Pa)')
    set(gca,'FontSize',18)
    figure(4)
    plot(TIME,stress(nb,:),'LineWidth',1)
    xlabel('TIME (h)')
    ylabel('Stress (Pa)')
    set(gca,'FontSize',18)
    figure(5)
    plot(TIME,stress(1,:),'LineWidth',1)
    xlabel('TIME (h)')
    ylabel('Stress (Pa)')
    set(gca,'FontSize',18)
    figure(6)
    plot(TIME,creepstrain(nb_c+1,:),'LineWidth',1)
    figure(7)
    plot(TIME,creepstrain(nb_c+nb_s+1,:),'LineWidth',1)
    figure(8)
    plot(TIME,creepstrain(nb,:),'LineWidth',1)
    
    
    
end
[m,n]=size(creepstrain);
fid=fopen('AlN_TZM_AlN_q10000_creepstrain_CTEofT.txt','wt');
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

[m,n]=size(stress);
fid=fopen('AlN_TZM_AlN_q10000_stress_CTEofT.txt','wt');
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

fid=fopen('AlN_TZM_AlN_q10000_Curvature_CTEofT.txt','wt');
fprintf(fid,'%12.5f\n',Curvature);
fclose(fid);

fid=fopen('AlN_TZM_AlN_q10000_TIMEN_CTEofT.txt','wt');
fprintf(fid,'%12.5f\n',TIME);
fclose(fid);
