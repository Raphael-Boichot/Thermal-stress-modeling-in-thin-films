function [vec_opt,fval]=GA(objectif,vec_ini,tol)
tic
% parametres
nb_pop = 10;           % nombre d'individus dans la population
nb_crois = 3;          % nombre des meilleurs individus autorisés à se croiser
p_mutation = 0.5;      % probabilité de mutation
p_croisement = 0.5;    % probabilité de croisement
borne = 0.1;           % intervalle de recherche x,y=[-val val]
max_steps = 15;        % nombre d'étapes maximales pour l'algo gen avant de relancer

options = optimset('display','off',...         % affichage des infos à chaque itération
    'largescale','on',...        % A préciser
    'MaxFunEvals',2000, ...      % nb max evaluation de fonction
    'MaxIter',    2000, ...      % nb max  d'itération
    'TolFun',     1e-12,  ...    % valeur critique fonction
    'TolX',       1e-12);        % tolérence itérés

%On récupère la taille du vecteur des données à optimiser
nb_var=vec_ini;
val=1;
% initialisation population et nb de cycle
vec_objectif=zeros(nb_pop,1);
pop=zeros(nb_pop, nb_var);

skip=0;
if isfile('sol.tmp')
    disp('Previous solution found, updating it')
    previous_sol=load('sol.tmp');
    [new_sol,fval]=fminsearch(objectif,previous_sol,options);
    if objectif(new_sol)<tol
        skip=1;
        disp('Skipping genetic algorithm, previous solutions was accurate enough after update')
        disp(['Obj_fun=',num2str(fval),' c=', num2str(new_sol(1,1)), ' tb=', num2str(new_sol(1,2)),' r=' , num2str(new_sol(1,3))]);
        vec_opt=new_sol;
    else
        disp('Previous solutions not accurate enough, starting from scratch')
    end
else
    disp('Previous solution not found, starting from scratch')
end

if skip==0;
    %création de la population intiale par algo de Monte Carlo
    disp('Starting with a Monte-Carlo algorithm')
    for i=1:1:20000
        for j=1:1:nb_var
            pop_ini(i,j)=borne*rand-2*borne*rand;
        end
        vecteur_initial(i,1)=objectif(pop_ini(i,:));
    end
    classement_ini=sortrows([pop_ini,vecteur_initial],nb_var+1);
    %On récupère les données classées
    pop=classement_ini(1:nb_pop,1:nb_var);

    disp('Starting genetic algorithm coupled to Nelder-Mead method')
    %toc
    obj=1e9;
    nb=1;
    % warning('off','all')
    while(tol<obj)
        nb=nb+1;
        % on evalue la population et on repère le meilleur et son indice
        % on stocke le meilleur dans un tableau pour tracer l'évolution
        % on créée une structure pour le classement des données
        for i=1:nb_pop
            %on balance un petit simplex au passage
            [pop(i,:)]=fminsearch(objectif,pop(i,:),options);
            vec_objectif(i) = objectif(pop(i,:));
        end

        classement=sortrows([pop,vec_objectif],nb_var+1);
        %On récupère les données classées
        pop=classement(:,1:nb_var);

        % on duplique le meilleur 1 fois dans la nouvelle population
        popnew(1,:)=pop(1,:);
        %popnew(2,:)=pop(1,:);

        %on stocke l'historique
        val(nb-1)=(objectif(pop(1,:)));

        % croisement
        for i=3:2:nb_pop
            % on selectionne 2 parents pour faire 2 enfants
            % choix parmis les 20 meilleurs
            i1=ceil(nb_crois*rand);
            i2=ceil(nb_crois*rand);
            %croisement et mutation
            popnew=croisement_GA(i,i1,i2,pop,popnew,p_croisement,p_mutation,borne);
        end
        pop=popnew;
        disp(['Epoch#',num2str(nb-1),' Obj_fun=',num2str(val(nb-1)),' c=', num2str(popnew(1,1)), ' tb=', num2str(popnew(1,2)),' r=' , num2str(popnew(1,3))]);
        obj=val(nb-1);
        if nb==max_steps
            obj=0;
            disp('Optimization failed, restarting from scratch');
        end
    end
    fval=obj;
    vec_opt=[popnew(1,:)];
end

save('sol.tmp','vec_opt','-ascii')
toc
