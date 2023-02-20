%Calculation of thermal stress in a multilayer coating
clc
clear all
close all
rng('shuffle', 'twister')
format 'ShortE'

global ilots;
global t_f;
global t_s;
global T_depot; 
global T_final;
global eps_thermique_film;
global eps_thermique_substrate;
global Y_substrate;
global Y_film;
global eps_misfit;
global eps_coalescence;

T_final=20;%°C température finale, peut être différente de T_initiale
T_depot=1200;
t_s=470e-6; %m, épaisseur du substrat
t_f=1e-6; %m, épaisseur du film
eps_misfit=-0.7e-2;
L_ilots=30e-9;
ilots=1;

%Data for AlN,
%alpha=(a0+a1(T°C)+a2(T°C)^2+a3(T°C)^3+a4(T°C)^4+a5(T°C)^5)*1e-6
E_film=340e9;      %Pa
nu_film=0.21;
a0_f=2.7405;
a1_f=1.1463e-2;
a2_f=-1.3902e-5;
a3_f=9.3822e-9;
a4_f=-3.3275e-12;
a5_f=4.8193e-16;
Y_film=E_film/(1-nu_film); %module biaxial en Pa
eps_coalescence=(2*L_ilots*(2*0.4-0)/Y_film)^0.5/L_ilots;

%Data for sapphire,
%alpha=(a0+a1(T°C)+a2(T°C)^2+a3(T°C)^3+a4(T°C)^4+a5(T°C)^5)*1e-6
E_substrate=431.24e9; %Pa
nu_substrate=0.285;
a0_s=4.844;
a1_s=1.4196e-2;
a2_s=-2.2829e-5;
a3_s=2.0237e-8;
a4_s=-8.5944e-12;
a5_s=1.3762e-15;
Y_substrate=E_substrate/(1-nu_substrate);

alpha_film = @(T)(a0_f+a1_f.*T+a2_f.*T.^2+a3_f.*T.^3+a4_f.*T.^4+a5_f.*T.^5)/1000000;
alpha_substrate = @(T)(a0_s+a1_s.*T+a2_s.*T.^2+a3_s.*T.^3+a4_s.*T.^4+a5_s.*T.^5)/1000000;

eps_thermique_film=quad(alpha_film,T_depot,T_final);
eps_thermique_substrate=quad(alpha_substrate,T_depot,T_final);

options = optimset('display','off',...         % affichage des infos à chaque itération
                   'largescale','off',...       % A préciser
                   'MaxFunEvals',2000, ...      % nb max evaluation de fonction
                   'MaxIter',    2000, ...      % nb max  d'itération
                   'TolFun',     1e-16,  ...    % valeur critique fonction
                   'TolX',       1e-16);         % tolérence itéréswarning('off','all')               
warning('off','all')

%Raw optimization with random search genetic algorithm
fval=1;
while fval>1e-8;
disp('Début de la recherche aléatoire');
f_ini=1e30;
for i=1:1:10000;
    x=[0.02.*rand-2.*0.02.*rand;-2.*rand.*t_s;10.*rand-2.*10.*rand];
    f=Contrainte_multicouche_Hsueh_BOICHOT(x);
    if f<f_ini;
        x_best=x;
        f_ini=f;
    end;
end;
disp('Fin de la recherche aléatoire');
x=x_best;

for i=1:1:4
[x,fval]=fminsearch('Contrainte_multicouche_Hsueh_BOICHOT',x,options);
end

if fval>1e-8;
disp('Echec de la méthode');
end


end
disp('Réussite de la méthode et convergence');
c=x(1);
t_b=x(2);
r=x(3);
disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);
if x(3)<0; disp('Dépôt en compression');end;
if x(3)>0; disp('Dépôt en tension');end;
disp(['Valeur de la fonction objectif=',num2str(fval)]);
epsilon=c+((t_f/2)-t_b)/r+eps_thermique_film-eps_misfit-eps_coalescence*ilots;
sigma_hsueh=-epsilon*Y_film;
disp(['Sigma méthode Hsueh=',num2str(sigma_hsueh, '%10.5e\n')]);
sigma_stoney=Y_substrate*t_s^2/(r*6*t_f);
disp(['Sigma approché formule Stoney=',num2str(sigma_stoney, '%10.5e\n')]);
raman_shift=657.67-sigma_hsueh*1e-9/4.04;
disp(['Raman shift Hsueh=',num2str(raman_shift)]);

