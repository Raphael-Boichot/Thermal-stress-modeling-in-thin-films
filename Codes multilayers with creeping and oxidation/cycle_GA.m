function[tao_accu,stress_cycle,creepstrain_cycle,Curvature_cycle]=cycle_GA(T_int,tao_int,stress_int,creepstrain_int);
global t_c;%thickness
global t_s;
global eps_thermique_coating_top;%misfit CTE
global eps_thermique_coating_bottom;
global eps_thermique_substrate;
global Y_substrate;
global Y_coating;
global nb_c;
global nb_s;
global nb;
global nb_t;
global deltat;
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
global T_elaboration;
global q_thermique;

%Data for substrate, TZM
E_substrate=320e9; %Pa
nu_substrate=0.3;
%nu_substrate=0;
Y_substrate=E_substrate/(1-nu_substrate);
%alpha_substrate=5.2e-6;%°C-1
a0_s=5.171085e+00;%a(TZM)=a(Mo), %polyfit order 5, Vasudevan1992
a1_s=-5.917661e-03;
a2_s=3.522412e-05;
a3_s=-5.012336e-08;
a4_s=3.164237e-11;
a5_s=-7.078432e-15;
lambda_substrate=103.5;%W/(m°C)
B_substrate=4.273e-12;%Pa-n.h-1
%B_substrate=0;
n_substrate=3.03;
R=8.314472;
Q_substrate=481038;%J/mol

%Data for coating, AlN
E_coating=340e9;%Pa
nu_coating=0.21;%http://accuratus.com/silicar.html
%nu_coating=0;
Y_coating=E_coating/(1-nu_coating);
%alpha_coating=4.5e-6;%°C-1
a0_c=2.47411e+00;
a1_c=1.11737e-02;
a2_c=-1.32814e-05;
a3_c=8.67958e-09;
a4_c=-2.94971e-12;
a5_c=4.05708e-16;
lambda_coating=30;%W/(m°C)
B_coating=exp(24.1);%MPa-ns-1
%B_coating=0;
n_coating=1;
Q_coating=586e3;

%thermal resistance
R_coating=t_c/lambda_coating;
R_substrate=t_s/lambda_substrate;
R_total=2*R_coating+R_substrate;

%input parameters
T_bottom=T_int(1);
T_ambiant=T_int(2);
T_top=T_bottom-R_total*q_thermique;

alpha_coating = @(T)(a0_c+a1_c.*T+a2_c.*T.^2+a3_c.*T.^3+a4_c.*T.^4+a5_c.*T.^5)/1000000;%polyfit by order 5
alpha_substrate = @(T)(a0_s+a1_s.*T+a2_s.*T.^2+a3_s.*T.^3+a4_s.*T.^4+a5_s.*T.^5)/1000000;%polyfit by order 5

z=zeros(nb,nb_t);
T=zeros(nb,nb_t);
stress_cycle=zeros(nb,nb_t);
creepstrain_cycle=zeros(nb,nb_t);
creepstrainrate=zeros(nb,(nb_t-1));

for j=1:nb_t%j:time
    if j==1
        T(:,j)=T_ambiant;
        stress_cycle(:,j)=stress_int;%initialization
        creepstrain_cycle(:,j)=creepstrain_int;%initialization

        for i=1:nb%z position
            if i<(nb_c+2)
                z(i)=t_s/2+t_c-t_c/nb_c*(i-1);
            end
            if (nb_c+1)<i && i<(nb_c+nb_s+2)
                z(i)=t_s/2-t_s/nb_s*(i-(nb_c+1));
            end
            if (nb_c+nb_s+1)<i
                z(i)=-t_s/2-t_c/nb_c*(i-(nb-nb_c));
            end
        end
    end
    if 1<j&&j<nb_t
        T(1,j)=T_top;
        T(nb,j)=T_bottom;
        time(j)=1/12+deltat*(j-2);%h
        %time(j)=1/12+deltat*sign(j-2)*(j-3)+1/1800*sign(j-2);
        tao(j-1)=deltat*(j-2)+tao_int;%holdtime accumulates
        %tao(j-1)=deltat*sign(j-2)*(j-3)+1/1800*sign(j-2);
    end

    if j==nb_t
        T(:,j)=T_ambiant;
        time(j)=time(j-1)+1/12;
        tao(j-1)=tao(j-2);
        tao_accu=tao(j-1);
    end
end


for j=2:nb_t
    z(:,j)=z(:,j-1);%no oxide growth during heating

    if j<nb_t%during heating
        for i=1:nb
            if i<(nb_c+1+1)%c2
                %T1(i,j)=T(1,j)-(T(1,j)-T(nb,j))/R_total*((t_c-z(i,j))/lambda_coating);
                T(i,j)=T(1,j)+q_thermique*((t_c+t_s/2-z(i,j))/lambda_coating);
                %creepstrain_cycle(i,j)=creepstrain_cycle(i,j-1)+creepstrainrate(i,j-1)*deltat;
            end
            if (nb_c+1)<i && i<(nb_c+nb_s+2)
                %T1(i,j)=T(1,j)-(T(1,j)-T(nb,j))/R_total*(-z(i,j)/lambda_substrate+R_coating);
                T(i,j)=T(1,j)+q_thermique*((t_s/2-z(i,j))/lambda_substrate+R_coating);
                %creepstrain_cycle(i,j)=creepstrain_cycle(i,j-1)+creepstrainrate(i,j-1)*deltat;
            end
            if (nb_c+nb_s+1)<i
                T(i,j)=T(1,j)+q_thermique*((-t_s/2-z(i,j))/lambda_coating+R_coating+R_substrate);
            end
        end
    end
    creepstrain_cycle(:,j)=creepstrain_cycle(:,j-1)+creepstrainrate(:,j-1)*deltat;
    %creepstrain_cycle(:,j)=creepstrain_cycle(:,j-1)+creepstrainrate(:,j-1)*(time(j)-time(j-1));%update creepstrain

    %if j==nb_t
    %growthstrainrate_thick(j)=0;
    %growthstrainrate_lat(j)=0;
    %    t_tgo(j)=t_tgo(j-1);
    %    z(i,j)=z(i,j-1);
    %end

    %R_tgo(j-1)=t_tgo(j)/lambda_tgo;
    %R_total(j-1)=R_coating_2+R_coating_1+R_tgo(j-1)+R_substrate;

    for i=1:nb
        if i<(nb_c+1+1)
            eps(i,j)=quad(alpha_coating,T_elaboration,T(i,j));
            %eps(i,j)=alpha_coating*(-T_elaboration+T(i,j));
        end
        if (nb_c+1)<i&&i<(nb_c+nb_s+2)
            eps(i,j)=quad(alpha_substrate,T_elaboration,T(i,j));
            %eps(i,j)=alpha_substrate*(-T_elaboration+T(i,j));
        end
        if (nb_s+nb_c+1)<i
            eps(i,j)=quad(alpha_coating,T_elaboration,T(i,j));
            %eps(i,j)=alpha_coating*(-T_elaboration+T(i,j));
        end
    end

    eps_thermique_coating_top=t_c/nb_c*sum(eps(2:(nb_c+1),j),1);
    %eps_thermique_coating_1=alpha_coating_1*t_c/nb_c*sum((T(2:(nb_c+1),j)-T_elaboration),1);%intégration alpha*deltaT*dz
    eps_thermique_substrate=t_s/nb_s*sum(eps((nb_c+1+1):(nb_c+nb_s+1),j),1);
    %eps_thermique_substrate_1=alpha_substrate_1*t_s/nb_s*sum((T((nb_c+1+1):nb,j)-T_elaboration),1);
    eps_thermique_coating_bottom=t_c/nb_c*sum(eps((nb_c+nb_s+1+1):nb,j),1);

    int_creep_substrate=sum(creepstrain_cycle((nb_c+1+1):(nb_c+nb_s+1),j),1)*t_s/nb_s;%intégration creepstrain*dz
    int_creep_coating_top=sum(creepstrain_cycle(2:(nb_c+1),j),1)*t_c/nb_c;
    int_creep_coating_bottom=sum(creepstrain_cycle((nb_c+nb_s+1+1):nb,j),1)*t_c/nb_c;

    moment_eps_coating_top=t_c/nb_c*sum(eps(2:(nb_c+1),j).*z(2:(nb_c+1),j),1);
    %moment_eps_coating_1=alpha_coating_1*t_c/nb_c*sum((T(2:(nb_c+1),j)-T_elaboration).*z(2:(nb_c+1),j),1);
    moment_eps_substrate=t_s/nb_s*sum(eps((nb_c+1+1):(nb_c+nb_s+1),j).*z((nb_c+1+1):(nb_c+nb_s+1),j),1);
    %moment_eps_substrate_1=alpha_substrate_1*t_s/nb_s*sum((T((nb_c+1+1):nb,j)-T_elaboration).*z((nb_c+1+1):nb,j),1);
    moment_eps_coating_bottom=t_c/nb_c*sum(eps((nb_c+nb_s+1+1):nb,j).*z((nb_c+nb_s+1+1):nb,j),1);

    moment_creep_coating_top=sum(creepstrain_cycle(2:(nb_c+1),j).*z(2:(nb_c+1),j),1)*t_c/nb_c;%intégration creepstrain*z*dz
    moment_creep_substrate=sum(creepstrain_cycle((nb_c+1+1):(nb_c+nb_s+1),j).*z((nb_c+1+1):(nb_c+nb_s+1),j),1)*t_s/nb_s;
    moment_creep_coating_bottom=sum(creepstrain_cycle((nb_c+nb_s+1+1):nb,j).*z((nb_c+nb_s+1+1):nb,j),1)*t_c/nb_c;

    fval=0;
    while fval==0
        [x,fval]=GA(@Contrainte_multicouche_Hsueh_4couches,3,1e-6);
    end

    c=x(1);
    t_b=x(2);
    r=x(3);

    Radis_cycle(j)=r;
    Curvature_cycle(j)=1/r;

    %if j<nb_t
    for i=1:nb
        if i<(nb_c+1+1)%c2
            stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN);
            %stress_cycle_1(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-alpha_coating_1*(T(i,j)-T_elaboration)-creepstrain_cycle(i,j)-strain_AlN);
            creepstrainrate(i,j)=B_coating*10^(-6*n_coating)*3600*(stress_cycle(i,j))^n_coating*exp(-Q_coating/(R*(273.15+T(i,j))));%1h=3600s
        end
        if (nb_c+1)<i&&i<(nb_c+nb_s+1+1)
            stress_cycle(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j));
            %stress_cycle_1(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-alpha_substrate_1*(T(i,j)-T_elaboration)-creepstrain_cycle(i,j));
            creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_substrate*(abs(stress_cycle(i,j)))^n_substrate*exp(-Q_substrate/(R*(273.15+T(i,j))));
            %steady state creep rate=B*stress^n*exp(-Q/(R*T)),T en K
        end
        if (nb_c+nb_s+1)<i
            stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN);
            creepstrainrate(i,j)=B_coating*10^(-6*n_coating)*3600*(stress_cycle(i,j))^n_coating*exp(-Q_coating/(R*(273.15+T(i,j))));
        end
    end
    creepstrainrate(:,nb_t-1)=0;
    creepstrainrate(:,nb_t)=0;
    %end

    %if j==nb_t
    %    for i=1:nb
    %        if i<(nb_c+1+1)%coating
    %            stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN);
    %            %stress_cycle_1(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-alpha_coating_1*(T(i,j)-T_elaboration)-creepstrain_cycle(i,j)-strain_AlN);
    %        else
    %        %substrate
    %            stress_cycle(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j));
    %            %stress_cycle_1(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-alpha_substrate*(T(i,j)-T_elaboration)-creepstrain_cycle(i,j));
    %        end
    %    end
    %end
    disp(['Step: ',num2str(j),'-----------------------------']);

end
disp('End of current temperature cycle');
%th_tgo_cycle=t_tgo;