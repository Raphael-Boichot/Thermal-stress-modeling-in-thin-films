%Calculation of thermal stress in a single layer coating on a substrate
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

T_final=20;%°C Final temperature, generally the room temperature
T_depot=1200;%Fabrication temperature
t_s=470e-6; %m, Substrate thickness
t_f=1e-6; %m, Film thickness
eps_misfit=-0.7e-2; %misfit stress
L_ilots=30e-9; %Size of initial nucleation islands in 3D growth
ilots=1; %0 to get rid of this effet

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
Y_film=E_film/(1-nu_film); %biaxial modulus Pa
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
Y_substrate=E_substrate/(1-nu_substrate);%biaxial modulus Pa

alpha_film = @(T)(a0_f+a1_f.*T+a2_f.*T.^2+a3_f.*T.^3+a4_f.*T.^4+a5_f.*T.^5)/1000000;
alpha_substrate = @(T)(a0_s+a1_s.*T+a2_s.*T.^2+a3_s.*T.^3+a4_s.*T.^4+a5_s.*T.^5)/1000000;

eps_thermique_film=quad(alpha_film,T_depot,T_final);
eps_thermique_substrate=quad(alpha_substrate,T_depot,T_final);

%Raw optimization with random search genetic algorithm
fval=0;
while fval==0;
[x,fval]=GA(@contrainte_multicouches,3,1e-8);
end;

disp('Global minimum found !');
c=x(1);%final stress
t_b=x(2);%neutral axis position (0 is the interface film/substrate)
r=x(3);%Radius of curvature
disp(['c=',num2str(c),'(-) tb=',num2str(t_b),'(m) r=',num2str(r),'(m)']);
if x(3)<0; disp('Compressive stress in "thin" film');end;
if x(3)>0; disp('Tensile stress in "thin" film');end;
disp(['Objective function at final step=',num2str(fval)]);
epsilon=c+((t_f/2)-t_b)/r+eps_thermique_film-eps_misfit-eps_coalescence*ilots;
sigma_hsueh=-epsilon*Y_film;
disp(['Sigma Hsueh algorithm=',num2str(sigma_hsueh, '%10.5e\n'), 'Pa']);
sigma_stoney=Y_substrate*t_s^2/(r*6*t_f);
disp(['Sigma with Stoney equation=',num2str(sigma_stoney, '%10.5e\n'), 'Pa']);
