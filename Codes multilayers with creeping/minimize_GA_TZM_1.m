%Calculation of thermal stress in a multilayer coating
function[r,stress_elaboration]=minimize_GA_TZM_1(T_elaboration,T_ambiant);


global t_c;
global t_s;
global T_depot; 
global T_final;
global strain_AlN;
global nb_c;
global nb_s;
global nb;


T_depot=T_elaboration;
T_final=T_ambiant;%°C température finale, peut être différente de T_initiale
%t_s=2000e-6; %m, épaisseur du substrat
%t_c=60e-6; %m, épaisseur du film

%Data for AlN,
%alpha=(a0+a1(T°C)+a2(T°C)^2+a3(T°C)^3+a4(T°C)^4+a5(T°C)^5)*1e-6
%E_film=340e9;      %Pa
%nu_film=0.21;
%a0_f=2.7405;
%a1_f=1.1463e-2;
%a2_f=-1.3902e-5;
%a3_f=9.3822e-9;
%a4_f=-3.3275e-12;
%a5_f=4.8193e-16;
%Y_film=E_film/(1-nu_film); %module biaxial en Pa


%Data for AlN, polycrystallin, 
%alpha=(a0+a1(T¡ãC)+a2(T¡ãC)^2+a3(T¡ãC)^3+a4(T¡ãC)^4+a5(T¡ãC)^5)*1e-6
E_film=340e9;      %Pa, http://accuratus.com/silicar.html
nu_film=0.21;      %http://accuratus.com/silicar.html
%nu_film=0;
%a0_f=2.47411e+00;
%a1_f=1.11737e-02;
%a2_f=-1.32814e-05;
%a3_f=8.67958e-09;
%a4_f=-2.94971e-12;
%a5_f=4.05708e-16;
Y_film=E_film/(1-nu_film); %module biaxial en Pa


%Data for TZM,
%alpha=(a0+a1(T°C)+a2(T°C)^2+a3(T°C)^3+a4(T°C)^4+a5(T°C)^5)*1e-6
E_substrate=320e9;
nu_substrate=0.3;
%nu_substrate=0;
%a0_s=5.171085e+00;%a(TZM)=a(Mo), %polyfit order 5, Vasudevan1992
%a1_s=-5.917661e-03;
%a2_s=3.522412e-05;
%a3_s=-5.012336e-08;
%a4_s=3.164237e-11;
%a5_s=-7.078432e-15;
Y_substrate=E_substrate/(1-nu_substrate);



%alpha_film = @(T)(a0_f+a1_f.*T+a2_f.*T.^2+a3_f.*T.^3+a4_f.*T.^4+a5_f.*T.^5)/1000000;
%alpha_substrate = @(T)(a0_s+a1_s.*T+a2_s.*T.^2+a3_s.*T.^3+a4_s.*T.^4+a5_s.*T.^5)/1000000;%polyfit by order 5

alpha_film=4.5e-6;
alpha_substrate=5.2e-6;

%eps_thermique_film=quad(alpha_film,T_depot,T_final);
%eps_thermique_substrate=quad(alpha_substrate,T_depot,T_final);

eps_thermique_film=alpha_film*(T_final-T_depot);
eps_thermique_substrate=alpha_substrate*(T_final-T_depot);


options = optimset('display','off',...         % affichage des infos ?chaque itération
                   'largescale','off',...       % A préciser
                   'MaxFunEvals',2000, ...      % nb max evaluation de fonction
                   'MaxIter',    2000, ...      % nb max  d'itération
                   'TolFun',     1e-16,  ...    % valeur critique fonction
                   'TolX',       1e-16);         % tolérence itérés               
warning('off','all')



%Raw optimization with random search genetic algorithm
%fval=0;
%while fval==0
%[x,fval]=GA_1(@Contrainte_multicouche_Hsueh_1couche,3,1e-8);
%end

r=inf;
t_b=0;
sum_biax_force=@(c) Y_substrate*t_s*(c-eps_thermique_substrate)+Y_film*t_c*(c-eps_thermique_film-strain_AlN)*2;

c=fzero(sum_biax_force,0);

disp('Réussite de la méthode et convergence');
%c=x(1);
%t_b=x(2);
%r=x(3);
disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);
%if x(3)>0 
%    disp('Dépôt en compression');
%end %definition de Hsueh: compression r>0
%if x(3)<0 
%    disp('Dépôt en tension');
%end

for i=1:nb%z position
	if i<(nb_c+2)
        z(i,1)=t_s/2+t_c-t_c/nb_c*(i-1);
        stress_elaboration(i,1)=Y_film*(c+(z(i,1)-t_b)/r-eps_thermique_film-strain_AlN);
    end
    if (nb_c+1)<i&&i<(nb_s+nb_c+2)
        z(i,1)=t_s/2-t_s/nb_s*(i-(nb_c+1));
        stress_elaboration(i,1)=Y_substrate*(c+(z(i,1)-t_b)/r-eps_thermique_substrate);
    end
    if (nb_s+nb_c+1)<i
        z(i,1)=-t_s/2-t_c/nb_c*(i-(nb-nb_c));
        stress_elaboration(i,1)=Y_film*(c+(z(i,1)-t_b)/r-eps_thermique_film-strain_AlN);
    end
end
disp('end')

%disp(['Valeur de la fonction objectif=',num2str(fval)]);
%epsilon=c+((t_f/2)-t_b)/r-eps_thermique_film-growthstrain;
%sigma_hsueh=epsilon*Y_film;
%disp(['Sigma m¨¦thode Hsueh=',num2str(sigma_hsueh, '%10.5e\n')]);
%epsilon1=c+(t_f-t_b)/r-eps_thermique_film-growthstrain;
%sigma_hsueh1=epsilon1*Y_film;
%sigma_stoney=Y_substrate*t_s^2/(-r*6*t_f);
%disp(['Sigma approch?formule Stoney=',num2str(sigma_stoney, '%10.5e\n')]);
%raman_shift=657.67-sigma_hsueh*1e-9*4.04;
%disp(['Raman shift Hsueh=',num2str(raman_shift)]);
%Y_substrate*(c+(-t_b)/r-eps_thermique_substrate)
%Y_substrate*(c+(-t_s-t_b)/r-eps_thermique_substrate)