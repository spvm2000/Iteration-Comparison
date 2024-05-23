function Data_o = WDO(Data_i,Data_o)
    t_c=clock;                             
    cnt=1;                               
    RT = 3;			
    g = 0.2;	    
    alp = 0.4;		
    c = 0.4;		
    maxV = 0.3.*Data_i.ub.*ones(1,Data_i.dim);	    
    minV = -0.3.*Data_i.ub.*ones(1,Data_i.dim);      
    IterCurve=zeros(1,Data_i.maxIter);
    pos = Data_i.X;
    V = initialization(Data_i.pop,maxV,minV,Data_i.dim);
    fitness = Data_i.F_value;

    [SortFitness,indexSort] = sort(fitness);
    gBest = pos(indexSort(1),:); 
    Best_Pos=indexSort(1);
    Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)=Best_Pos;
    gBestFitness = SortFitness(1);  

    pos = pos(indexSort,:);
    fitness = SortFitness;

    for t = 1:Data_i.maxIter     
       for i = 1:Data_i.pop

            a = randperm(Data_i.dim);    
            velot = V(i,a);

            V(i,:) = (1- alp)*V(i,:)-(g*pos(i,:))+ ...
				    abs(1-1/i)*((gBest-pos(i,:)).*RT)+ ...
				    (c*velot/i);
 
            V(i,:) = BoundaryCheck(V(i,:),maxV,minV,Data_i.dim);
  
            pos(i,:) = pos(i,:) + V(i,:);

            pos(i,:) = BoundaryCheck(pos(i,:),Data_i.ub,Data_i.lb,Data_i.dim);
     
            fitness(i) = Data_i.fobj(pos(i,:));
            if fitness(i)<gBestFitness
                gBestFitness = fitness(i);
                gBest = pos(i,:);
            end 
       end
       [SortFitness,indexSort] = sort(fitness);
       Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)=gBestFitness;
    
       pos = pos(indexSort,:);   
       V = V(indexSort,:);
       fitness = SortFitness;
       Data_o.Best_Pos(Data_i.now_fun_index,Data_i.now_test_iter)= indexSort(1);             
       IterCurve(t)=gBestFitness;
end

        Data_o.exe_time(Data_i.now_fun_index,Data_i.now_test_iter)=etime(clock,t_c);              
        Data_o.min_t(Data_i.now_fun_index,Data_i.now_test_iter)=t;                            
    Data_o.Best_fitness(Data_i.now_fun_index,Data_i.now_test_iter)= gBestFitness;

Data_o.IterCurve(Data_i.now_fun_index,Data_i.now_test_iter,:)=IterCurve;
end