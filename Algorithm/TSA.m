function Data_o = TSA(Data_i,Data_o)
    t_c=clock;                             
    cnt=1;                              
    low = ceil(0.1*Data_i.pop);
    high = ceil(0.25*Data_i.pop);
    ST = 0.1;
    trees = Data_i.X;
    fitness = Data_i.F_value;
    IterCurve=zeros(1,Data_i.maxIter);
    [SortFitness,indexSort] = sort(fitness);
    Best_Pos=indexSort(1);
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Pos;
    gBest = trees(indexSort(1),:);  
    gBestFitness = SortFitness(1);  
    for t = 1:Data_i.maxIter
       for i = 1:Data_i.pop
            seedNum =fix(low+(high-low)*rand)+1;   
            seeds=zeros(seedNum,Data_i.dim);
            obj_seeds=zeros(1,seedNum);
            [minimum,min_indis]=min(fitness); 
            bestParams=trees(min_indis,:);
            for j = 1:seedNum
                komsu=fix(rand*Data_i.pop)+1;
                while(i==komsu) 
                   komsu=fix(rand*Data_i.pop)+1;
                end
                seeds(j,:)=trees(j,:);
                for d=1:Data_i.dim
                    if(rand<ST)
                        seeds(j,d)=trees(i,d)+(bestParams(d)-trees(komsu,d))*(rand-0.5)*2;
                    else
                        seeds(j,d)=trees(i,d)+(trees(i,d)-trees(komsu,d))*(rand-0.5)*2;
                    end
                end
                seeds(j,:) = BoundaryCheck(seeds(j,:),Data_i.ub,Data_i.lb,Data_i.dim);
                obj_seeds(j)=Data_i.fobj(seeds(j,:));
            end
            [mintohum,mintohum_indis]=min(obj_seeds);
            if(mintohum<fitness(i))
                trees(i,:)=seeds(mintohum_indis,:);
                fitness(i)=mintohum;
            end     
       end
       [min_tree,min_tree_index]=min(fitness);
       if(min_tree<gBestFitness)
            gBestFitness=min_tree;
            gBest=trees(min_tree_index,:);
            Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=min_tree_index;
       end
       Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=gBestFitness;
       IterCurve(t)=gBestFitness;
    end
Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                             
Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end