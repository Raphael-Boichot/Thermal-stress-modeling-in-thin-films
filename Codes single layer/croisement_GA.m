function popnew=croisement_alternate(i,i1,i2,pop,popnew,p_croisement,p_mutation,borne)
     
    [test,size_vec]=size(pop(1,:));

    chromosome_p1='';
    chromosome_p2='';
    chiffre_hex='0123456789abcdef';
    p_mutation=p_mutation*rand;
    p_croisement=p_croisement*rand;
     %On créé un chromosome à partir des gènes des deux parents
     for k=1:1:size_vec;
            chromosome_p1=[chromosome_p1,num2hex(pop(i1,k))];
            chromosome_p2=[chromosome_p2,num2hex(pop(i2,k))];
     end;
     
     chromosome_e1=chromosome_p1;
     chromosome_e2=chromosome_p2;
          
     %On croise les chromosomes avec une probabilité p_croisement
     sens_lecture=1;
     
     [test,longueur]=size(chromosome_p1);
     
     for j=1:1:longueur;
         
     if p_croisement>rand; 
     sens_lecture= sens_lecture*-1;
     end;
         
     if sens_lecture==1;
     chromosome_e1(j)=chromosome_p1(j);
     chromosome_e2(j)=chromosome_p2(j);
     end;
     
     if sens_lecture==-1;   
     chromosome_e1(j)=chromosome_p2(j);
     chromosome_e2(j)=chromosome_p1(j);
     end
     
     %si on rencontre un chiffre, on mute les deux chromosomes
     if p_mutation>rand;
   
     mutation1=chiffre_hex(ceil(16*rand));
     mutation2=chiffre_hex(ceil(16*rand)); 
     chromosome_e1(j)=mutation1;
     chromosome_e2(j)=mutation2;

     end
     end

     %On vérifie la viabilité des chromosomes
     
     for k=1:1:size_vec;
         
         pos=(k-1)*16+1;
         vec_sortie1(k)=hex2num(chromosome_e1(pos:pos+15));         
         vec_sortie2(k)=hex2num(chromosome_e2(pos:pos+15)); 
         
         if abs(vec_sortie1(k)>1e15); vec_sortie1(k)=borne*rand-2*borne*rand;end;
         if abs(vec_sortie2(k)>1e15); vec_sortie2(k)=borne*rand-2*borne*rand;end;
         
         if abs(vec_sortie1(k)<1e-15); vec_sortie1(k)=borne*rand-2*borne*rand;end;
         if abs(vec_sortie2(k)<1e-15); vec_sortie2(k)=borne*rand-2*borne*rand;end;

     end;

     popnew(i,:)=vec_sortie1;
     popnew(i+1,:)=vec_sortie2;

     
    
        
     
     
     
     