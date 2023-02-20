function [vec_opt,fval]=GA(objectif,vec_ini,tol);
tic
% parametres
 nb_pop = 20;           % nb dans la population
 nb_crois = 5;          % nb meilleurs autoris�s � se croiser
 p_mutation = 0.5;      % proba de mutation
 p_croisement = 0.5;    % proba de croisement
 borne = 0.1;           % intervalle de recherche x,y=[-val val]

 
%On r�cup�re la taille du vecteur des donn�es � optimiser
nb_var=vec_ini;
val=1;
% initialisation population et nb de cycle
vec_objectif=zeros(nb_pop,1);
%val=zeros(nb_cycle,1);
pop=zeros(nb_pop, nb_var);

options = optimset('display','off',...         % affichage des infos � chaque it�ration
                   'largescale','off',...       % A pr�ciser
                   'MaxFunEvals',2000, ...      % nb max evaluation de fonction
                   'MaxIter',    2000, ...      % nb max  d'it�ration
                   'TolFun',     1e-16,  ...    % valeur critique fonction
                   'TolX',       1e-16);         % tol�rence it�r�s  

%cr�ation de la population intiale
for i=1:1:nb_pop;
    for j=1:1:nb_var;
    pop(i,j)=borne*rand-2*borne*rand;
    end;
end;
obj=1e9;
nb=1;
    while(tol<obj)
nb=nb+1;
    % on evalue la population et on rep�re le meilleur et son indice
    % on stocke le meilleur dans un tableau pour tracer l'�volution
    %On cr��e une structure pour le classement des donn�es
    for i=1:1:nb_pop;
        %on balance un petit simplex au passage
        [pop(i,:)]=fminsearch(objectif,pop(i,:),options);
        vec_objectif(i) = objectif(pop(i,:));
    end;
    classement=sortrows([pop,vec_objectif],nb_var+1);
    %On r�cup�re les donn�es class�es
    pop=classement(:,[1:nb_var]);

    % on duplique le meilleur 2 fois dans la nouvelle population
    %remplissage ligne 1
    popnew(1,:)=pop(1,:);
    popnew(2,:)=pop(1,:);
    
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

    if rem(nb,1)==0;
    disp(['Iteration #',num2str(nb-1),' Obj function=',num2str(val(nb-1))]);
    disp(popnew(1,:));
    end;
    obj=val(nb-1);
    if nb==10;
        obj=0;
        disp('Algorithm failed to find global minimum, reload');
    end;
    end
    
    % on affiche les r�sultats  
    
    vec_opt=[popnew(1,:)];
    fval=obj;
toc