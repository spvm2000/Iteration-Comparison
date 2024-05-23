function Data_o = SOS(Data_i,Data_o)
rand_num=[];                         
ti=clock;                               
eco = Data_i.X;                        
fitness=Data_i.F_value;              
[Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter),Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)]=min(Data_i.F_value);
IterCurve=zeros(1,Data_i.maxIter);
n=Data_i.dim;       ecosize=Data_i.pop;     maxFE=Data_i.maxIter;  fobj=Data_i.fobj;
FE=1;       ub=Data_i.ub;       lb=Data_i.lb;
while FE<=maxFE 
    for i=1:ecosize
        [bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);
        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        mutualVector=mean([eco(i,:);eco(j,:)]);
        BF1=round(1+rand); BF2=round(1+rand);
        ecoNew1=eco(i,:)+rand(1,n).*(bestOrganism-BF1.*mutualVector); 
        ecoNew2=eco(j,:)+rand(1,n).*(bestOrganism-BF2.*mutualVector);
        ecoNew1=BoundaryCheck(ecoNew1,ub,lb,n); 
        ecoNew2=BoundaryCheck(ecoNew2,ub,lb,n);
        fitnessNew1=fobj(ecoNew1);
        fitnessNew2=fobj(ecoNew2);

        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end
        if fitnessNew2<fitness(j)
           fitness(j)=fitnessNew2;
           eco(j,:)=ecoNew2;
        end


        j=i;
        while i==j
            seed=randperm(ecosize); 
            j=seed(1);                  
        end
        ecoNew1=eco(i,:)+(rand(1,n)*2-1).*(bestOrganism-eco(j,:));
        ecoNew1=BoundaryCheck(ecoNew1,ub,lb,n);
        fitnessNew1=fobj(ecoNew1);
        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end

        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        parasiteVector=eco(i,:);
        seed=randperm(n);           
        pick=seed(1:ceil(rand*n));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        fitnessParasite=fobj(parasiteVector);

        if fitnessParasite < fitness(j)
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
        end
    end
    [value,index]=min(fitness);
    if Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)>value
        Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=value;
        Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=index;
    end    
    IterCurve(FE)=value;
    FE=FE+1;
end    
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,ti);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=FE;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end