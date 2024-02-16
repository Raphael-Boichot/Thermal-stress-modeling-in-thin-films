function[T,T_oxi_out,z_out,tao_accu,th_cycle,growthstrain_lat_top_accu,growthstrain_lat_bottom_accu,eps_tgo_top_out,eps_tgo_bottom_out,stress_cycle,creepstrain_cycle,Curvature_cycle]=cycle_GA(T_oxi_int,z_int,tao_int,th_int,growthstrain_lat_top_int,growthstrain_lat_bottom_int,eps_tgo_top_int,eps_tgo_bottom_int,stress_int,creepstrain_int)

global t_s;%thickness
global t_preoxide;
global t_int;
global t_tgo_top;
global t_tgo_bottom;
global t_c_top;
global t_c_bottom;
global Y_substrate;
global Y_coating;
global Y_tgo;
global Y_oxide;
global alpha_coating;
global alpha_substrate;
global alpha_tgo;
global alpha_oxide;
global nb_c1;
global nb_c2;
global nb_preoxide
global nb_s;
global nb_o;
global nb;
global nb_t;
global growthstrain_lat_top;
global growthstrain_lat_bottom;
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
global j;
global k;
global deltat;
global strain_AlN;
global T_elaboration;
global T_preoxidation;
global T_ambiant;
global size_z_sub;
global size_z_preoxide;
global size_z_coating;
global size_z_tgo;
global C;
global Lateral_preoxide;

%Data for substrate, APMT
lambda_substrate=27;%W/(m°C)
B_substrate=78.978;%MPa-n.s-1
n_substrate=4.9827;
R=8.314472;
Q_substrate=453e3;%J/mol

%Data for coating, AlN
lambda_coating=30;%W/(m°C)
B_coating=exp(24.1)*200;%MPa-n.s-1
n_coating=1;
Q_coating=586e3;
%Q_coating=625e3;
M_coating=40.9882;%g/mol
density_coating=3.255;%g/cm3
k_0=35.9202068;%mg/cm2/h
Q_oxidation=11705*R;%J/mol
delta_M_oxidation=3/2*16*2-14*2;

%Data for tgo, Al2O3
lambda_tgo=10;%W/(m°C)
B_tgo=exp(15.68301);%MPa-n.s-1
n_tgo=1.64812;
Q_tgo=460e3;
deltaEd=230.89e3;%J/mol
thetat=5.7;%s
density_tgo=3.97;%g/cm3
M_tgo=101.9613;%g/mol
%B_tgo=3e-9;%MPa-n.s-1
%B_tgo=3e-10;
%n_tgo=1;
%Q_tgo=0;

%Data for preoxide,Al2O3
lambda_preoxide=10;%W/(m°C)
B_preoxide=exp(15.68301);%MPa-n.s-1
n_preoxide=1.64812;
Q_preoxide=460e3;
%B_preoxide=3e-9;%MPa-n.s-1
%B_preoxide=3e-10;%MPa-n.s-1
%n_preoxide=1;
%Q_preoxide=0;

%thermal resistance
R_substrate=t_s/lambda_substrate;
R_preoxide=t_preoxide/lambda_preoxide;

z=zeros(nb,nb_t);
T=zeros(nb,nb_t);
stress_cycle=zeros(nb,nb_t);
creepstrain_cycle=zeros(nb,nb_t);
creepstrainrate=zeros(nb,(nb_t-1));

for j=1:nb_t%j:time
    %initialization
    if j==1
        T(:,j)=T_ambiant;
        z(:,j)=z_int;
        t_tgo_top(j)=th_int(1,1);
        t_tgo_bottom(j)=th_int(2,1);
        t_c_top(j)=th_int(3,1);
        t_c_bottom(j)=th_int(4,1);
        t_tgo_top1=t_int(1);
        t_tgo_bottom1=t_int(2);
        t_c_top1=t_int(3);
        t_c_bottom1=t_int(4);
        growthstrain_lat_top(j)=growthstrain_lat_top_int;
        growthstrain_lat_bottom(j)=growthstrain_lat_bottom_int;
        stress_cycle(:,j)=stress_int;
        creepstrain_cycle(:,j)=creepstrain_int;

        R_coating_top(j)=t_c_top(j)/lambda_coating;
        R_coating_bottom(j)=t_c_bottom(j)/lambda_coating;
        R_tgo_top(j)=t_tgo_top(j)/lambda_tgo;
        R_tgo_bottom(j)=t_tgo_bottom(j)/lambda_tgo;
        R_total(j)=R_coating_top(j)+R_coating_bottom(j)+R_tgo_top(j)+R_tgo_bottom(j)+R_substrate+2*R_preoxide;
    end

    if 1<j&&j<nb_t
        %T(nb,j)=T_bottom;
        time(j)=1/12+deltat*(j-2);%h
        tao(j-1)=deltat*(j-2)+tao_int;%holdtime accumulates
    end

    if j==nb_t
        T(:,j)=T_ambiant;
        time(j)=time(j-1)+1/12;
        tao(j-1)=tao(j-2);
        tao_accu=tao(j-1);
    end
end

%no oxide growth during heating
t_tgo_top(2)=t_tgo_top(1);
t_tgo_bottom(2)=t_tgo_bottom(1);
t_c_top(2)=t_c_top(1);
t_c_bottom(2)=t_c_bottom(1);
z(:,2)=z(:,1);
growthstrain_lat_top(2)=growthstrain_lat_top(1);
growthstrain_lat_bottom(2)=growthstrain_lat_bottom(1);
R_coating_top(2)=R_coating_top(1);
R_coating_bottom(2)=R_coating_bottom(1);
R_tgo_top(2)=R_tgo_top(1);
R_tgo_bottom(2)=R_tgo_bottom(1);
R_total(2)=R_total(1);

nb_tgo_top(2)=ceil(t_tgo_top(2)/size_z_tgo);
nb_tgo_bottom(2)=ceil(t_tgo_bottom(2)/size_z_tgo);
nb_c_top(2)=ceil(t_c_top(2)/size_z_coating);
nb_c_bottom(2)=ceil(t_c_bottom(2)/size_z_coating);


for j=2:nb_t

    if j<nb_t
        %update Temperature
        for i=1:nb
            %    if i<=(nb_o+1)%tgo_top
            %        if t_tgo_top(j)==0%T(i,j)=T(1,j) if no oxide scale;
            %            T(i,j)=T_bottom-R_total(j)*q_thermique;
            %        else
            %thermal conduction
            %            T(i,j)=T(nb,j)-q_thermique*((z(i,j)-t_s/2-t_preoxide-t_c_top(j))/lambda_tgo+...
            %                   R_coating_top(j)+R_substrate+2*R_preoxide+R_coating_bottom(j)+R_tgo_bottom(j));
            %        end
            %    end
            %    if (nb_o+1)<i&&i<=(nb_c1+nb_o+1)%c_top
            %        T(i,j)=T(nb,j)-q_thermique*((z(i,j)-t_s/2-t_preoxide)/lambda_coating+R_substrate+...
            %               2*R_preoxide+R_coating_bottom(j)+R_tgo_bottom(j));
            %    end
            %    if (nb_o+nb_c1+1)<i&&i<=(nb_c1+nb_o+nb_preoxide+1)%preoxide_top
            %        T(i,j)=T(nb,j)-q_thermique*((z(i,j)-t_s/2)/lambda_preoxide+R_substrate+R_preoxide+...
            %               R_coating_bottom(j)+R_tgo_bottom(j));
            %    end
            %    if (nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+nb_preoxide+nb_c1+nb_o+1)%substrate
            %        T(i,j)=T(nb,j)-q_thermique*((z(i,j)-(-t_s/2))/lambda_substrate+R_preoxide+...
            %               R_coating_bottom(j)+R_tgo_bottom(j));
            %    end
            %    if (nb_s+nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+2*nb_preoxide+nb_c1+nb_o+1)%preoxide_bottom
            %        T(i,j)=T(nb,j)-q_thermique*((z(i,j)-(-t_s/2-t_preoxide))/lambda_preoxide+R_coating_bottom(j)+...
            %               R_tgo_bottom(j));
            %    end
            %    if (nb_s+2*nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+2*nb_preoxide+nb_c1+nb_c2+nb_o+1)%c_bottom
            %        T(i,j)=T(nb,j)-q_thermique*((z(i,j)-(-t_s/2-t_preoxide-t_c_bottom(j)))/lambda_coating+R_tgo_bottom(j));
            %    end
            %    if (nb-nb_o)<i%tgo_bottom
            %        if t_tgo_bottom(j)==0
            %            T(i,j)=T(nb,j);%T(i,j)=T(nb,j) if no oxide scale;
            %        else
            %            T(i,j)=T(nb,j)-q_thermique*(z(i,j)-(-t_s/2-t_preoxide-t_c_bottom(j)-t_tgo_bottom(j)))/lambda_tgo;
            %        end
            %    end
            T(:,j)=1100;
        end
    end

    %no oxide scale growth during cooling
    z(:,nb_t)=z(:,(nb_t-1));

    %update creepstrain
    creepstrain_cycle(:,j)=creepstrain_cycle(:,j-1)+creepstrainrate(:,j-1)*deltat;%update creepstrain

    %store oxidation temperature for c_top
    if k==1
        T_oxi_top=T(nb_o+1,2);
    else
        T_oxi_top=T_oxi_int;
    end
    T_oxi_bottom=T(nb,2);

    if j<nb_t
        %update oxide scale thickness
        %avancement_tgo_top=(tao(j)*3600*1/thetat*exp(-deltaEd/(R*(T_oxi_top+273.15))))^0.5;
        %avancement_tgo_bottom=(tao(j)*3600*1/thetat*exp(-deltaEd/(R*(T(nb,j)+273.15))))^0.5;
        %t_tgo_top(j+1)=t_tgo_top1+avancement_tgo_top*t_c_top1*(density_coating*M_tgo)/(density_tgo*M_coating)/4;
        %t_tgo_bottom(j+1)=t_tgo_bottom1+avancement_tgo_bottom*t_c_bottom1*(density_coating*M_tgo)/(density_tgo*M_coating)/4;
        %t_c_top(j+1)=t_c_top1*(1-0.5*avancement_tgo_top);
        %t_c_bottom(j+1)=(1-0.5*avancement_tgo_bottom)*t_c_bottom1;
        t_tgo_top(j+1)=1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T_oxi_top+273.15)))*M_tgo/density_tgo/delta_M_oxidation;
        t_tgo_bottom(j+1)=1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T(nb,j)+273.15)))*M_tgo/density_tgo/delta_M_oxidation;
        t_c_top(j+1)=t_c_top1-1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T_oxi_top+273.15)))*M_tgo/density_tgo/delta_M_oxidation*(2*M_coating/density_coating/(M_tgo/density_tgo))^(1/3);
        t_c_bottom(j+1)=t_c_bottom1-1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T(nb,j)+273.15)))*M_tgo/density_tgo/delta_M_oxidation*(2*M_coating/density_coating/(M_tgo/density_tgo))^(1/3);
        %t_c_top(j+1)=t_c_top1-1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T_oxi_top+273.15)))*2*M_coating/density_coating/delta_M_oxidation;
        %t_c_bottom(j+1)=t_c_bottom1-1e-5*tao(j)*k_0*exp(-Q_oxidation/(R*(T(nb,j)+273.15)))*2*M_coating/density_coating/delta_M_oxidation;

        %growthstrain_lat(j+1)=0.0019*tao(j)^0.53655;
        if j==(nb_t-1)
            growthstrain_lat_top(j+1)=growthstrain_lat_top(j);
            growthstrain_lat_bottom(j+1)=growthstrain_lat_bottom(j);
        else
            %growthstrain_lat_top_rate=1/2*lateral_strain_constant_AlN*...
            %                          (1/thetat*exp(-deltaEd/(R*(T_oxi_top+273.15))))^0.5*...
            %                          t_c_top1*(density_coating*M_tgo)/(density_tgo*M_coating)/4*...
            %                          (tao(j)*3600)^(-1/2);%s-1
            %growthstrain_lat_bottom_rate=1/2*lateral_strain_constant_AlN*...
            %                             (1/thetat*exp(-deltaEd/(R*(T(nb,j)+273.15))))^0.5*...
            %                             t_c_bottom1*(density_coating*M_tgo)/(density_tgo*M_coating)/4*...
            %                             (tao(j)*3600)^(-1/2);%s-1
            %growthstrain_lat_top_rate=3600*1e-5*lateral_strain_constant_AlN*k_0*exp(-Q_oxidation/(R*T_oxi_top))*...
            %                          M_tgo/density_tgo/delta_M_oxidation;
            %growthstrain_lat_bottom_rate=3600*1e-5*lateral_strain_constant_AlN*k_0*exp(-Q_oxidation/(R*T(nb,j)))*...
            %                             M_tgo/density_tgo/delta_M_oxidation;
            %growthstrain_lat_top(j+1)=growthstrain_lat_top(j)+growthstrain_lat_top_rate*tao(j)*3600;
            %growthstrain_lat_bottom(j+1)=growthstrain_lat_bottom(j)+growthstrain_lat_bottom_rate*tao(j)*3600;
            growthstrain_lat_top(j+1)=1e-3;
            growthstrain_lat_bottom(j+1)=1e-3;
        end

        %update thermal resistance
        R_coating_top(j+1)=t_c_top(j+1)/lambda_coating;
        R_coating_bottom(j+1)=t_c_bottom(j+1)/lambda_coating;
        R_tgo_top(j+1)=t_tgo_top(j+1)/lambda_tgo;
        R_tgo_bottom(j+1)=t_tgo_bottom(j+1)/lambda_tgo;
        R_total(j+1)=R_coating_top(j+1)+R_coating_bottom(j+1)+R_tgo_top(j+1)+R_tgo_bottom(j+1)+R_substrate+2*R_preoxide;

        %update z
        nb_tgo_top(j+1)=ceil(t_tgo_top(j+1)/size_z_tgo);
        nb_tgo_bottom(j+1)=ceil(t_tgo_bottom(j+1)/size_z_tgo);
        nb_c_top(j+1)=ceil(t_c_top(j+1)/size_z_coating);
        nb_c_bottom(j+1)=ceil(t_c_bottom(j+1)/size_z_coating);
        for i=1:nb
            if i<=(nb_o+1)
                %T1(i,j)=T(i,j);
                if i<=(nb_o-nb_tgo_top(j+1)+1)
                    z(i,j+1)=t_tgo_top(j+1)+t_c_top(j+1)+t_s/2+t_preoxide;
                else
                    z(i,j+1)=t_c_top(j+1)+t_preoxide+t_s/2+size_z_tgo*(nb_o+1-i);
                end
            end
            if (nb_o+1)<i&&i<=(nb_c1+nb_o+1)%c_top
                if i<=(nb_c1+nb_o+1-nb_c_top(j+1))
                    z(i,j+1)=t_preoxide+t_s/2+t_c_top(j+1);
                else
                    z(i,j+1)=t_preoxide+t_s/2+size_z_coating*(nb_o+nb_c1+1-i);
                end
                %T1(i,j)=T(nb,j)-q_thermique*((z(i,j+1)-t_s/2)/lambda_coating+R_substrate+R_coating_bottom(j)+R_tgo_bottom(j));
                %delta z change for t+deltat, temperature at t for point at z(t+deltat)
            end
            if (nb_o+nb_c1+1)<i&&i<=(nb_c1+nb_o+nb_preoxide+1)%preoxide_top
                z(i,j+1)=z(i,j);
            end
            if (nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+nb_preoxide+nb_c1+nb_o+1)%substrate
                z(i,j+1)=z(i,j);
                %T1(i,j)=T(1,j)+q_thermique*((t_s/2-z(i,j+1))/lambda_substrate+R_coating_top(j)+R_tgo_top(j));
                %T1(i,j)=T(nb,j)-q_thermique*((z(i,j)-(-t_s/2))/lambda_substrate+R_coating_bottom(j)+R_tgo_bottom(j));
                %TT(i-(nb_c+nb_o+1),j)=T1(i,j)-T(i,j);%for test, TT=0 for every time
            end
            if (nb_s+nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+2*nb_preoxide+nb_c1+nb_o+1)%preoxide_bottom
                z(i,j+1)=z(i,j);
            end
            if (nb_s+2*nb_preoxide+nb_c1+nb_o+1)<i&&i<=(nb_s+2*nb_preoxide+nb_c1+nb_c2+nb_o+1)%c_bottom
                if i<=(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j+1))
                    z(i,j+1)=-t_s/2-t_preoxide-size_z_coating*(i-(nb_o+nb_c1+2*nb_preoxide+nb_s+1));
                else
                    z(i,j+1)=-t_s/2-t_preoxide-t_c_bottom(j+1);
                end
                %T1(i,j)=T(nb,j)-q_thermique*((z(i,j+1)-(-t_s/2-t_c_bottom(j)))/lambda_coating+R_tgo_bottom(j));
                %delta z change for t+deltat, temperature at t for point at z(t+deltat)
            end
            if (nb-nb_o)<i%tgo_bottom
                %T1(i,j)=T(i,j);
                if i>=(nb_o+nb_c1+nb_c2+2*nb_preoxide+nb_s+nb_tgo_bottom(j+1)+1)
                    z(i,j+1)=-t_tgo_bottom(j+1)-t_s/2-t_preoxide-t_c_bottom(j+1);
                else
                    z(i,j+1)=-t_s/2-t_preoxide-t_c_bottom(j+1)-size_z_tgo*(i-(nb-nb_o));
                end
            end
        end
        %i_tgo_top(j)=ceil((t_tgo_top(j+1)-t_tgo_top(j))/(t_tgo_top(j+1)/nb_o));
        %i_tgo_bottom(j)=nb-ceil((t_tgo_bottom(j+1)-t_tgo_bottom(j))/(t_tgo_bottom(j+1)/nb_o))+1;
    end

    %i_tgo_top(1)=i_tgo_top(2);
    %i_tgo_bottom(1)=i_tgo_bottom(2);
    if k==1
        nb_tgo_top(2)=0;
        nb_tgo_bottom(2)=0;
        nb_c_top(2)=nb_c1;
        nb_c_bottom(2)=nb_c2;
    end


    if j==nb_t
        %t_tgo_top(j)=t_tgo_top(j-1);
        %t_tgo_bottom(j)=t_tgo_bottom(j-1);
        %t_c_top(j)=t_c_top(j-1);
        %t_c_bottom(j)=t_c_bottom(j-1);
        growthstrain_lat_top_accu=growthstrain_lat_top(j);
        growthstrain_lat_bottom_accu=growthstrain_lat_bottom(j);
        %T1(:,j)=T(:,j);
        %i_tgo_top(j)= i_tgo_top(j-1);
        %i_tgo_bottom(j)= i_tgo_bottom(j-1);
        %z(:,j)=z(:,(j-1));
    end

    %thermal expansion calculation
    for i=1:nb
        if i<(nb_o+2)%tgo_top
            eps(i,j)=integral(alpha_tgo,T_oxi_top,T(i,j));
            %eps1(i,j)=integral(alpha_tgo,T_oxi_top,T1(i,j));
        end
        if (nb_o+1)<i&&i<=(nb_c1+nb_o+1)%c_top
            eps(i,j)=integral(alpha_coating,T_elaboration,T(i,j));
            %eps1(i,j)=integral(alpha_coating,T_elaboration,T1(i,j));%eps at t for point at z(t+deltat)
        end
        if (nb_o+nb_c1+1)<i&&i<=(nb_c1+nb_o+nb_preoxide+1)%preoxide_top
            eps(i,j)=integral(alpha_oxide,T_preoxidation,T(i,j));
        end
        if (nb_preoxide+nb_c1+nb_o+1)<i&&i<(nb_s+nb_preoxide+nb_c1+nb_o+1)%substrate
            eps(i,j)=integral(alpha_substrate,T_preoxidation,T(i,j));
            %eps1(i,j)=integral(alpha_substrate,T_elaboration,T1(i,j));
        end
        if (nb_s+nb_preoxide+nb_c1+nb_o+1)<=i&&i<(nb_s+2*nb_preoxide+nb_c1+nb_o+1)%preoxide_bottom
            eps(i,j)=integral(alpha_oxide,T_preoxidation,T(i,j));
            %eps1(i,j)=integral(alpha_coating,T_elaboration,T1(i,j));%eps at t for point at z(t+deltat)
        end
        if (nb_s+2*nb_preoxide+nb_c1+nb_o+1)<=i&&i<(nb_s+2*nb_preoxide+nb_c1+nb_c2+nb_o+1)%c_bottom
            eps(i,j)=integral(alpha_coating,T_elaboration,T(i,j));
            %eps1(i,j)=integral(alpha_coating,T_elaboration,T1(i,j));%eps at t for point at z(t+deltat)
        end
        if (nb-nb_o)<=i
            eps(i,j)=integral(alpha_tgo,T_oxi_bottom,T(i,j));
            %eps1(i,j)=integral(alpha_tgo,T_bottom,T1(i,j));
        end
    end

    %at the beginning of the first cycle, no oxide (k=1,j=2)
    if eps_tgo_top_int==0&&eps_tgo_bottom_int==0
        %    eps_tgo_top1=integral(alpha_coating,T_elaboration,T_top);
        %    eps_tgo_bottom1=eps((nb-nb_o),2);
        eps_tgo_top1=eps_tgo_top_int;
        eps_tgo_bottom1=eps_tgo_bottom_int;
        %else
        %    eps_tgo_top1=eps_tgo_top_int;
        %    eps_tgo_bottom1=eps_tgo_bottom_int;
    end

    eps_thermique_tgo_top=size_z_tgo*sum(eps((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j),1)+...
        eps((nb_o+1-nb_tgo_top(j)+1),j)*(z((nb_o+1-nb_tgo_top(j)),j)-z((nb_o+1-nb_tgo_top(j)+1),j));

    eps_thermique_coating_top=size_z_coating*sum(eps((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j),1)+...
        eps((nb_o+nb_c1+1-nb_c_top(j)+1),j)*(z((nb_o+nb_c1+1-nb_c_top(j)),j)-z((nb_o+nb_c1+1-nb_c_top(j)+1),j));

    eps_thermique_preoxide_top=size_z_preoxide*sum(eps((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j),1)+...
        eps((nb_o+nb_c1+1+1),j)*(t_preoxide+t_s/2-z((nb_o+nb_c1+1+1),j));

    eps_thermique_substrate=size_z_sub*sum(eps((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j),1)+...
        eps((nb_preoxide+nb_c1+nb_o+1+1),j)*(t_s/2-z((nb_preoxide+nb_c1+nb_o+1+1),j))+...
        eps((nb_s+nb_preoxide+nb_c1+nb_o),j)*(t_s/2+z((nb_s+nb_preoxide+nb_c1+nb_o),j));

    eps_thermique_preoxide_bottom=size_z_preoxide*sum(eps((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j),1)+...
        eps((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    eps_thermique_preoxide_bottom1=size_z_preoxide*sum(eps((nb_s+nb_preoxide+nb_c1+nb_o+1+1):(nb_s+2*nb_preoxide+nb_c1+nb_o),j),1)+...
        eps((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    eps_thermique_coating_bottom=size_z_coating*sum(eps((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j),1)+...
        eps((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*(z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)-z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)+1),j));

    eps_thermique_tgo_bottom=size_z_tgo*sum(eps((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j),1)+...
        eps(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*(z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)-z((nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)+1),j));

    moment_eps_tgo_top=size_z_tgo*sum(eps((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j).*z((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j),1)+...
        eps((nb_o+1-nb_tgo_top(j)+1),j)*z((nb_o+1-nb_tgo_top(j)+1),j)*(z((nb_o+1-nb_tgo_top(j)),j)-z((nb_o+1-nb_tgo_top(j)+1),j));

    moment_eps_coating_top=size_z_coating*sum(eps((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j).*z((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j),1)+...
        eps((nb_o+nb_c1+1-nb_c_top(j)+1),j)*z((nb_o+nb_c1+1-nb_c_top(j)+1),j)*(z((nb_o+nb_c1+1-nb_c_top(j)),j)-z((nb_o+nb_c1+1-nb_c_top(j)+1),j));

    moment_eps_preoxide_top=size_z_preoxide*sum(eps((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j).*z((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j),1)+...
        eps((nb_o+nb_c1+1+1),j)*z((nb_o+nb_c1+1+1),j)*(t_preoxide+t_s/2-z((nb_o+nb_c1+1+1),j));

    moment_eps_substrate=size_z_sub*sum(eps((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j).*z((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j),1)+...;
        eps((nb_preoxide+nb_c1+nb_o+1+1),j)*z((nb_preoxide+nb_c1+nb_o+1+1),j)*(t_s/2-z((nb_preoxide+nb_c1+nb_o+1+1),j))+...
        eps((nb_s+nb_preoxide+nb_c1+nb_o),j)*z((nb_s+nb_preoxide+nb_c1+nb_o),j)*(t_s/2+z((nb_s+nb_preoxide+nb_c1+nb_o),j));

    moment_eps_preoxide_bottom=size_z_preoxide*sum(eps((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j).*z((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j),1)+...
        eps((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*z((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    moment_eps_coating_bottom=size_z_coating*sum(eps((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j).*z((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j),1)+...
        eps((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*(z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)-z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)+1),j));

    moment_eps_tgo_bottom=size_z_tgo*sum(eps((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j).*z((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j),1)+...
        eps(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*(z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)-z((nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)+1),j));

    int_creep_tgo_top=size_z_tgo*sum(creepstrain_cycle((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j),1)+...
        creepstrain_cycle((nb_o+1-nb_tgo_top(j)+1),j)*(z((nb_o+1-nb_tgo_top(j)),j)-z((nb_o+1-nb_tgo_top(j)+1),j));

    int_creep_coating_top=size_z_coating*sum(creepstrain_cycle((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+1-nb_c_top(j)+1),j)*(z((nb_o+nb_c1+1-nb_c_top(j)),j)-z((nb_o+nb_c1+1-nb_c_top(j)+1),j));

    int_creep_preoxide_top=size_z_preoxide*sum(creepstrain_cycle((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+1+1),j)*(t_preoxide+t_s/2-z((nb_o+nb_c1+1+1),j));

    int_creep_substrate=size_z_sub*sum(creepstrain_cycle((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j),1)+...
        creepstrain_cycle((nb_preoxide+nb_c1+nb_o+1+1),j)*(t_s/2-z((nb_preoxide+nb_c1+nb_o+1+1),j))+...
        creepstrain_cycle((nb_s+nb_preoxide+nb_c1+nb_o),j)*(t_s/2+z((nb_s+nb_preoxide+nb_c1+nb_o),j));

    int_creep_preoxide_bottom=size_z_preoxide*sum(creepstrain_cycle((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j),1)+...
        creepstrain_cycle((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    int_creep_preoxide_bottom1=size_z_preoxide*sum(creepstrain_cycle((nb_s+nb_preoxide+nb_c1+nb_o+1+1):(nb_s+2*nb_preoxide+nb_c1+nb_o),j),1)+...
        creepstrain_cycle((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    int_creep_coating_bottom=size_z_coating*sum(creepstrain_cycle((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*(z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)-z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)+1),j));

    int_creep_tgo_bottom=size_z_tgo*sum(creepstrain_cycle((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j),1)+...
        creepstrain_cycle(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*(z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)-z((nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)+1),j));

    moment_creep_tgo_top=size_z_tgo*sum(creepstrain_cycle((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j).*z((nb_o+1-nb_tgo_top(j)+1+1):(nb_o+1),j),1)+...
        creepstrain_cycle((nb_o+1-nb_tgo_top(j)+1),j)*z((nb_o+1-nb_tgo_top(j)+1),j)*(z((nb_o+1-nb_tgo_top(j)),j)-z((nb_o+1-nb_tgo_top(j)+1),j));

    moment_creep_coating_top=size_z_coating*sum(creepstrain_cycle((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j).*z((nb_o+nb_c1+1-nb_c_top(j)+1+1):(nb_o+nb_c1+1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+1-nb_c_top(j)+1),j)*z((nb_o+nb_c1+1-nb_c_top(j)+1),j)*(z((nb_o+nb_c1+1-nb_c_top(j)),j)-z((nb_o+nb_c1+1-nb_c_top(j)+1),j));

    moment_creep_preoxide_top=size_z_preoxide*sum(creepstrain_cycle((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j).*z((nb_o+nb_c1+1+1+1):(nb_c1+nb_o+nb_preoxide+1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+1+1),j)*z((nb_o+nb_c1+1+1),j)*(t_preoxide+t_s/2-z((nb_o+nb_c1+1+1),j));

    moment_creep_substrate=size_z_sub*sum(creepstrain_cycle((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j).*z((nb_preoxide+nb_c1+nb_o+1+1+1):(nb_s+nb_preoxide+nb_c1+nb_o),j),1)+...;
        creepstrain_cycle((nb_preoxide+nb_c1+nb_o+1+1),j)*z((nb_preoxide+nb_c1+nb_o+1+1),j)*(t_s/2-z((nb_preoxide+nb_c1+nb_o+1+1),j))+...
        creepstrain_cycle((nb_s+nb_preoxide+nb_c1+nb_o),j)*z((nb_s+nb_preoxide+nb_c1+nb_o),j)*(t_s/2+z((nb_s+nb_preoxide+nb_c1+nb_o),j));

    moment_creep_preoxide_bottom=size_z_preoxide*sum(creepstrain_cycle((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j).*z((nb_s+nb_preoxide+nb_c1+nb_o+1):(nb_s+2*nb_preoxide+nb_c1+nb_o-1),j),1)+...
        creepstrain_cycle((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*z((nb_s+2*nb_preoxide+nb_c1+nb_o),j)*(t_s/2+t_preoxide+z((nb_s+2*nb_preoxide+nb_c1+nb_o),j));

    moment_creep_coating_bottom=size_z_coating*sum(creepstrain_cycle((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j).*z((nb_s+2*nb_preoxide+nb_c1+nb_o+1):(nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)-1),j),1)+...
        creepstrain_cycle((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)*(z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)),j)-z((nb_o+nb_c1+2*nb_preoxide+nb_s+nb_c_bottom(j)+1),j));

    moment_creep_tgo_bottom=size_z_tgo*sum(creepstrain_cycle((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j).*z((nb-nb_o):(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)-1),j),1)+...
        creepstrain_cycle(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)*(z(nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j),j)-z((nb_o+nb_preoxide*2+nb_c1+nb_c2+nb_s+nb_tgo_bottom(j)+1),j));

    if (t_c_top(j)-t_c_bottom(j))==0
        t_b=0;
        r=Inf;
        %sum of biaxial stress=0
        sum_biaxial_stress = @(c) Y_tgo*(c*t_tgo_top(j)-eps_tgo_top1*t_tgo_top(j)-eps_thermique_tgo_top-int_creep_tgo_top-t_tgo_top(j)*growthstrain_lat_top(j))+...
            Y_coating*(c*t_c_top(j)-eps_thermique_coating_top-int_creep_coating_top-t_c_top(j)*strain_AlN-t_c_top(j)*C)+...
            Y_oxide*(c*t_preoxide-eps_thermique_preoxide_top-int_creep_preoxide_top-t_preoxide*Lateral_preoxide)+...
            Y_substrate*(t_s*c-eps_thermique_substrate-int_creep_substrate)+...
            Y_oxide*(c*t_preoxide-eps_thermique_preoxide_bottom-int_creep_preoxide_bottom-t_preoxide*Lateral_preoxide)+...
            Y_coating*(c*t_c_bottom(j)-eps_thermique_coating_bottom-int_creep_coating_bottom-t_c_bottom(j)*strain_AlN-t_c_bottom(j)*C)+...
            Y_tgo*(c*t_tgo_bottom(j)-eps_tgo_bottom1*t_tgo_bottom(j)-eps_thermique_tgo_bottom-int_creep_tgo_bottom-t_tgo_bottom(j)*growthstrain_lat_bottom(j));
        c=fzero(sum_biaxial_stress,0);
        disp(['c=',num2str(c),' tb=',num2str(t_b),' r=',num2str(r)]);
    else
        fval=0;
        while fval==0;
            [x,fval]=GA(@Contrainte_multicouche_Hsueh_cycle,3,1e-6);
        end
        disp('Converging, next step beginning');
        c=x(1);
        t_b=x(2);
        r=x(3);
    end

    if j==2
        if eps_tgo_top_int==0&&eps_tgo_bottom_int==0
            eps_tgo_top1=c+(z(nb_o+1,j)-t_b)/r;
            eps_tgo_bottom1=c+(z(nb-nb_o,j)-t_b)/r;
            eps_tgo_top_int=eps_tgo_top1;
            eps_tgo_bottom_int=eps_tgo_bottom1;
        else
            eps_tgo_top1=eps_tgo_top_int;
            eps_tgo_bottom1=eps_tgo_bottom_int;
        end
    end

    Radis_cycle(j)=r;
    Curvature_cycle(j)=1/r;

    if j<nb_t
        for i=1:nb
            if i<=(nb_o+1)%tgo_top
                if t_tgo_top(j)==0
                    stress_cycle(i,j)=0;
                    creepstrainrate(i,j)=0;
                else
                    if i<(nb_o-nb_tgo_top(j)+1)
                        stress_cycle(i,j)=0;
                        creepstrainrate(i,j)=0;
                    else
                        stress_cycle(i,j)=Y_tgo*(c+(z(i,j)-t_b)/r-eps_tgo_top1-eps(i,j)-creepstrain_cycle(i,j)-growthstrain_lat_top(j));
                        creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_tgo*10^(-6*n_tgo)*3600*(abs(stress_cycle(i,j)))^n_tgo*exp(-Q_tgo/(R*(273.15+T(i,j))));
                    end
                end
            end

            if (nb_o+1)<i&&i<=(nb_o+nb_c1+1)%ctop
                stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN-C);
                creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_coating*10^(-6*n_coating)*3600*(abs(stress_cycle(i,j)))^n_coating*exp(-Q_coating/(R*(273.15+T(i,j))));%1h=3600s
            end

            if (nb_o+nb_c1+1)<i&&i<=(nb_o+nb_c1+nb_preoxide+1)%preoxide_top
                stress_cycle(i,j)=Y_oxide*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-Lateral_preoxide);
                creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_preoxide*10^(-6*n_preoxide)*3600*(abs(stress_cycle(i,j)))^n_preoxide*exp(-Q_preoxide/(R*(273.15+T(i,j))));
            end
            if (nb_o+nb_c1+nb_preoxide+1)<i&&i<(nb_o+nb_c1+nb_s+nb_preoxide+1)%substrate
                stress_cycle(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j));
                creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_substrate*10^(-6*n_substrate)*3600*(abs(stress_cycle(i,j)))^n_substrate*exp(-Q_substrate/(R*(273.15+T(i,j))));
                %steady state creep rate=B*stress^n*exp(-Q/(R*T)),T en K
            end
            if (nb_o+nb_c1+nb_s+nb_preoxide+1)<=i&&i<(nb_o+nb_c1+nb_s+2*nb_preoxide+1)%preoxide_bottom
                stress_cycle(i,j)=Y_oxide*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-Lateral_preoxide);
                creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_preoxide*10^(-6*n_preoxide)*3600*(abs(stress_cycle(i,j)))^n_preoxide*exp(-Q_preoxide/(R*(273.15+T(i,j))));
            end
            if (2*nb_preoxide+nb_s+nb_c1+nb_o+1)<=i&&i<(nb-nb_o)%c_bottom
                stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN-C);
                creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_coating*10^(-6*n_coating)*3600*(abs(stress_cycle(i,j)))^n_coating*exp(-Q_coating/(R*(273.15+T(i,j))));
            end
            if (nb-nb_o)<=i
                if t_tgo_bottom(j)==0
                    stress_cycle(i,j)=0;
                    creepstrainrate(i,j)=0;
                else
                    if i>(nb_tgo_bottom(j)+nb_o+nb_c1+nb_c2+nb_s+nb_preoxide*2+1)
                        stress_cycle(i,j)=0;
                        creepstrainrate(i,j)=0;
                    else
                        stress_cycle(i,j)=Y_tgo*(c+(z(i,j)-t_b)/r-eps_tgo_bottom1-eps(i,j)-creepstrain_cycle(i,j)-growthstrain_lat_bottom(j));
                        creepstrainrate(i,j)=sign(stress_cycle(i,j))*B_tgo*10^(-6*n_tgo)*3600*(abs(stress_cycle(i,j)))^n_tgo*exp(-Q_tgo/(R*(273.15+T(i,j))));
                    end
                end
            end
        end
        creepstrainrate(:,nb_t-1)=0;
    end

    disp(['----------------step time ',num2str(j),'----------------'])

    if j==(nb_t-1)
        disp('End of holding');
    end

    if j==nb_t
        for i=1:nb
            if i<=(nb_o+1)%tgo_top
                if t_tgo_top(j)==0
                    stress_cycle(i,j)=0;
                else
                    if i<(nb_o-nb_tgo_top(j)+1)
                        stress_cycle(i,j)=0;
                    else
                        stress_cycle(i,j)=Y_tgo*(c+(z(i,j)-t_b)/r-eps_tgo_top1-eps(i,j)-creepstrain_cycle(i,j)-growthstrain_lat_top(j));
                    end
                end
            end

            if (nb_o+1)<i&&i<=(nb_o+nb_c1+1)%ctop
                stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN-C);
            end

            if (nb_o+nb_c1+1)<i&&i<=(nb_o+nb_c1+nb_preoxide+1)%preoxide_top
                stress_cycle(i,j)=Y_oxide*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-Lateral_preoxide);
            end
            if (nb_o+nb_c1+nb_preoxide+1)<i&&i<(nb_o+nb_c1+nb_s+nb_preoxide+1)%substrate
                stress_cycle(i,j)=Y_substrate*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j));
            end
            if (nb_o+nb_c1+nb_s+nb_preoxide+1)<=i&&i<(nb_o+nb_c1+nb_s+2*nb_preoxide+1)%preoxide_bottom
                stress_cycle(i,j)=Y_oxide*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-Lateral_preoxide);
            end
            if (2*nb_preoxide+nb_s+nb_c1+nb_o+1)<=i&&i<(nb-nb_o)%c_bottom
                stress_cycle(i,j)=Y_coating*(c+(z(i,j)-t_b)/r-eps(i,j)-creepstrain_cycle(i,j)-strain_AlN-C);
            end
            if (nb-nb_o)<=i
                if t_tgo_bottom(j)==0
                    stress_cycle(i,j)=0;
                else
                    if i>(nb_tgo_bottom(j)+nb_o+nb_c1+nb_c2+nb_s+nb_preoxide*2+1)
                        stress_cycle(i,j)=0;
                    else
                        stress_cycle(i,j)=Y_tgo*(c+(z(i,j)-t_b)/r-eps_tgo_bottom1-eps(i,j)-creepstrain_cycle(i,j)-growthstrain_lat_bottom(j));
                    end
                end
            end
        end
    end
end
disp('End of one thermal cycle');
T_oxi_out=T_oxi_top;
th_cycle(1,:)=t_tgo_top;
th_cycle(2,:)=t_tgo_bottom;
th_cycle(3,:)=t_c_top;
th_cycle(4,:)=t_c_bottom;
z_out=z(:,nb_t);
eps_tgo_top_out=eps_tgo_top1;
eps_tgo_bottom_out=eps_tgo_bottom1;
